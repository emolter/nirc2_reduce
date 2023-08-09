import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits
from .image import Image
from .prettycolors import make_colormap, get_colormap
from scipy.signal import medfilt
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import rotate
import astroscrappy
from image_registration.chi2_shifts import chi2_shift
from image_registration.fft_tools.shift import shiftnd, shift2d
import importlib.resources
import warnings
import yaml
import nirc2_reduce.data.header_kw_dicts as inst_info


'''
To do
-----
weight stacks by integration time (?)
'''


def crop_center(frame, subc):
    """
    Parameters
    ----------
    frame : np.array, required.
    subc : int, required.
    
    Returns
    -------
    np.array
    """
    shp = np.min([frame.shape[0], frame.shape[1]])
    ctr = int(shp / 2)
    ll = int(ctr - subc / 2)
    ul = int(ctr + subc / 2)
    frame = frame[ll:ul, ll:ul]
    return frame


class Observation:
    """
    Description
    -----------
    Generic class containing dithers and nods
    To define a custom dither pattern, inherit from this class
    See Bxy3() and Nod() objects for usage
    """

    def __init__(self, fnames, instrument):
        """
        Parameters
        ----------
        fnames : list or str, required. 
            filenames of data frames. 
            if str, assume passing single data frame
        instrument : str, required.
            name of instrument used, e.g. nirc2. 
            will look for file data/{instrument}.yaml
            in order to scrub header keywords
        
        Attributes
        ----------
        dummy_fits : nirc2_reduce.image.Image of the zeroth frame
        frames : np.array, shape (n, x, y) for the n input fits images
        instrument : input instrument string
        subc : int. size of the subarray
        target : str. the object you observed, scrubbed from the header.
        """
        self.instrument = instrument.lower()
        with importlib.resources.open_binary(inst_info, f"{self.instrument}.yaml") as file:
            yaml_bytes = file.read()
            self.header_kw_dict = yaml.safe_load(yaml_bytes)
        
        if type(fnames) == str:
            fnames = [
                fnames,
            ]  # needed to make pass through all the for loops

        self.dummy_fits = Image(fnames[0])  # used to hijack header info
        self.frames = np.asarray([Image(f).data for f in fnames])
        self.subc = int(self.dummy_fits.header[self.header_kw_dict['subc']])
        self.target = self.dummy_fits.header[self.header_kw_dict['object']]

        if self.frames[0].shape[0] != self.subc or self.frames[0].shape[1] != self.subc:
            # subarrays smaller than 512x512 on Keck are not square. chop out center
            frames_square = []
            for frame in self.frames:
                frames_square.append(crop_center(frame, self.subc))
            self.frames = np.asarray(frames_square)

    def apply_flat(self, fname):
        """
        Description
        -----------
        applies flatfield correction
        
        Parameters
        ----------
        fname : str, required. 
            filename of flat
        """
        flat = Image(fname).data
        # handle subarrays
        if flat.shape[0] != self.subc:
            flat = crop_center(flat, self.subc)

        frames_flat = []
        for frame in self.frames:
            with np.errstate(divide="ignore", invalid="ignore"):
                frames_flat.append(frame / flat)
        self.frames = np.asarray(frames_flat)

    def apply_badpx_map(self, fname, kernel_size=7):
        """
        Description
        -----------
        Replaces all bad pixels in input map with the median of the pixels
        around it.
        
        Parameters
        ----------
        fname : str, required.
            filename of bad pixel fits file to apply
        kernel_size : int, optional.
            how large a region around which to take median
        """
        badpx_map = Image(fname).data
        # handle subarrays
        if badpx_map.shape[0] != self.subc:
            badpx_map = crop_center(badpx_map, self.subc)
        bad_indices = np.where(badpx_map == 0)
        ## test case - kernel of 7 shown to remove all bad pixels including weird stripe in bottom left of detector
        # smoothed_badpx_map = medfilt(badpx_map,kernel_size = 7)
        # fig, (ax0,ax1) = plt.subplots(1,2, figsize=(10,6))
        # ax0.imshow(badpx_map, origin = 'lower')
        # ax1.imshow(smoothed_badpx_map, origin = 'lower')
        # plt.show()
        frames_badpx = []
        for frame in self.frames:
            smoothed = medfilt(frame, kernel_size=kernel_size)
            frame[bad_indices] = smoothed[bad_indices]
            frames_badpx.append(frame)
        self.frames = np.asarray(frames_badpx)

    def dewarp(self):
        """
        Description
        -----------
        Apply nirc2 distortion correction using updated maps 
        from Ghez, Lu galactic center group (Service et al. 2016)
        doi:10.1088/1538-3873/128/967/095004
        Spline interpolation used to interpolate original image. 
        The spline is then evaluated at new pixel locations computed from the map.
        Custom distortion solutions should be placed in the 
        nirc2_reduce/nirc2_reduce/data/ folder
        
        Parameters
        ----------
        warpx_file : str, optional.
            distortion map to apply in x
        warpy_file : str, optional.
            distortion map to apply in y
        
        Notes
        -----
        From Service+16 "These are lookup tables generated by evaluating the fits 
        from the previous section at the center of every pixel on the 
        NIRC2 detector and are the values that should be added to 
        raw NIRC2 positions to shift them to a distortion-free frame"
        So the position offsets in the fits files should be added. 
        That this works right was checked by looking at difference maps 
        in a subdirectory of tests/
        and comparing with the expected vectors in Service+16
        """
        warpx_file = self.header_kw_dict['warpx_file']
        warpy_file = self.header_kw_dict['warpy_file']
        distortion_source_x = importlib.resources.open_binary(
            "nirc2_reduce.data", f"{warpx_file}"
        )
        distortion_source_y = importlib.resources.open_binary(
            "nirc2_reduce.data", f"{warpy_file}"
        )
        warpx = fits.getdata(distortion_source_x)
        warpy = fits.getdata(distortion_source_y)

        # handle subarrays
        if (
            self.subc != warpx.shape[0]
        ):  # hardcode because dewarp arrays will always be full size of detector
            ctr = warpx.shape[0]/2
            ll = int(ctr - self.subc / 2)
            ul = int(ctr + self.subc / 2)
            warpx = warpx[ll:ul, ll:ul]
            warpy = warpy[ll:ul, ll:ul]
        # plt.imshow(warpx,origin='lower')
        # plt.show()
        szx = self.frames[0].shape[0]
        szy = self.frames[0].shape[0]
        xx = np.linspace(0, szx - 1, szx)
        yy = np.linspace(0, szy - 1, szy)
        x, y = np.meshgrid(xx, yy)
        mapx = x - warpy  # check if plus or minus!!
        mapy = y + warpx  # check if plus or minus!!
        flatx = mapx.flatten()
        flaty = mapy.flatten()
        frames_dewarp = []
        for frame in self.frames:
            frame_spline = RectBivariateSpline(xx, yy, frame)
            frame_dw = frame_spline.ev(flatx, flaty).reshape(frame.shape).T
            frames_dewarp.append(frame_dw)
            # plt.imshow(frame_dw, origin='lower')
            # plt.show()
        self.frames = np.asarray(frames_dewarp)

    def rotate(self, custom_angle=0.0, beta=0.252):
        """
        Description
        -----------
        Rotate frame to North up plus an additional custom angle
        if header contains 'ROTPOSN' and 'INSTANGL' keywords (i.e., nirc2-like),
        then rotation is custom_angle + (ROTPOSN-INSTANGL) - beta, counterclockwise.
        otherwise, rotation is simply custom_angle - beta, counterclockwise.
        see https://github.com/jluastro/nirc2_distortion/wiki
        
        Parameters
        ----------
        custom_angle : float, optional. default 0.0. units degrees.
            user-defined additional angle to rotate frames
            counterclockwise is positive to agree with scipy.ndimage.rotate
        beta : float, optional. default 0.252 (nirc2 post 2015). units degrees.
            rotation required to get North up according to astrometry solution
            clockwise is positive to agree with definition of beta given by
            the nirc2 distortion wiki
            
        Notes
        -----
        This should be generalized to include scopes other than nirc2 more easily
        however, for typical observations ROTPOSN - INSTANGL is zero anyway
        """
        hdr = self.dummy_fits.header
        try:
            total_rotation_ccw = custom_angle + hdr["ROTPOSN"] - hdr["INSTANGL"] - beta
        except KeyError:
            warnings.warn(
                "header keywords ROTPOSN and/or INSTANGL not found; setting rotation to custom_angle - beta"
            )
        self.frames = np.array(
            [rotate(frame, total_rotation_ccw, reshape=False) for frame in self.frames]
        )

    def remove_cosmic_rays(self):
        """
        Description
        -----------
        Detects cosmic rays using the astroscrappy package
        """
        frames_cosray = []
        for frame in self.frames:
            crmask, cleanarr = astroscrappy.detect_cosmics(frame, cleantype="medmask")
            # fig, (ax0,ax1,ax2) = plt.subplots(1,3, figsize=(15,6))
            # ax0.imshow(frame,origin = 'lower')
            # ax1.imshow(crmask,origin = 'lower')
            # ax2.imshow(cleanarr,origin = 'lower')
            # plt.show()
            frames_cosray.append(cleanarr)
        self.frames = np.asarray(frames_cosray)

    def per_second(self):
        """
        Description
        -----------
        Changes units to counts/second by dividing by
        ITIME x COADDS
        
        Parameters
        ----------
        itime_kw : str, optional.
            keyword to scrub fits header for integration time
        coadd_kw : str, optional.
            keyword to scrub fits header for coadds
        """
        header = self.dummy_fits.header
        tint = header[self.header_kw_dict['itime']]
        coadd = header[self.header_kw_dict['coadds']]
        persec_frames = []
        for frame in self.frames:
            persec_frames.append(frame / (tint * coadd))
        self.frames = np.asarray(persec_frames)

    def apply_photometry_frames(self, flux_per):
        """
        Description
        -----------
        Simply multiplies each frame by flux density / (cts/s)
        
        Parameters
        ----------
        flux_per : float, required. units [flux density / (cts/s)]
             conversion factor
        """
        self.frames = np.asarray([frame * flux_per for frame in self.frames])

    def apply_photometry_final(self, flux_per):
        """
        Description
        -----------
        Simply multiplies each frame by flux density / (cts/s)
        
        Parameters
        ----------
        flux_per : float, required. units [flux density / (cts/s)]
            conversion factor 
        """
        self.final = self.final * flux_per

    def stack(self):
        """
        Description
        -----------
        Cross-correlate the images applying sub-pixel shift.
        Shift found using DFT upsampling method from image_registration.
        Stack them on top of each other to increase SNR.
        """
        shifted_data = [self.frames[0]]
        for frame in self.frames[1:]:
            [dx, dy, dxerr, dyerr] = chi2_shift(self.frames[0], frame)
            # error is nonzero only if you include per-pixel error of each image as an input. Should eventually do that, but no need for now.
            shifted = shift2d(frame, -1 * dx, -1 * dy)
            shifted_data.append(shifted)

        self.final = np.median(shifted_data, axis=0)

    def crop_final(self, bw):
        """
        Description
        -----------
        Applies crop to the final image. 
        
        Parameters
        ----------
        bw : int, required. 
            border width
        """
        szx, szy = self.final.shape[0], self.final.shape[1]
        self.final = self.final[bw : szx - bw, bw : szy - bw]

    def plot_frames(self, png_file=None):
        """
        Description
        -----------
        Plot the individual frames any step in the process
        
        Parameters
        ----------
        png_file : str or None, optional.
            filename to which image is saved
            if None, image not saved
        """
        n = len(self.frames)
        fig, axes = plt.subplots(1, n, figsize=(3 * n, 4))
        for i in range(len(axes)):
            plotframe = self.frames[i]
            vmax = np.nanmax(medfilt(plotframe, kernel_size=7))
            ax = axes[i]
            ax.imshow(plotframe, origin="lower", vmin=0, vmax=vmax)
            ax.set_title("Frame %d" % i)
            ax.set_xticks([])
            ax.set_yticks([])
        if png_file is not None:
            fig.savefig(png_file, dpi=300)
        plt.show()

    def plot_final(self, show=True, png_file=None):
        """
        Parameters
        ----------
        show : bool, optional.
            show plot?
        png_file : str or None, optional.
            filename to which image is saved
            if None, image not saved
        """
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        try:
            cmap = get_colormap(self.target.split(" ")[0])
        except:
            print("No custom colormap defined for target, setting to default")
            cmap = cm.viridis
        ax.imshow(self.final, cmap=cmap, origin="lower")
        plt.tight_layout()
        if png_file is not None:
            fig.savefig(png_file, dpi=300)
        if show:
            plt.show()
        plt.close()

    def write_frames(self, outfiles):
        """
        Description
        -----------
        Write each individual frame to .fits at any step in the process
        
        Parameters
        ----------
        outfiles : list, required. 
            filenames to write
        """
        for i in range(len(self.frames)):
            frame = self.frames[i]
            outfile = outfiles[i]
            hdulist_out = self.dummy_fits.hdulist
            hdulist_out[0].data = frame
            hdulist_out[0].writeto(outfile, overwrite=True)

    def write_final(self, outfile):
        """
        Description
        -----------
        Writes the final stacked frame to .fits
        
        Parameters
        ----------
        outfile : str, required. 
            filename to write
        png_file : str or None, optional.
            filename to which image is saved
            if None, image not saved
        """
        hdulist_out = self.dummy_fits.hdulist
        hdulist_out[0].header["OBJECT"] = self.target + "_STACKED"
        hdulist_out[0].data = self.final
        hdulist_out[0].writeto(outfile, overwrite=True)


class Nod(Observation):
    """
    Description
    -----------
    simple on-off nod dither pattern
    
    To do
    -----
    how to make it so self.final gets defined if there is no stack() function?
    """

    def __init__(self, data_fname, sky_fname, instrument):
        """
        Parameters
        ----------
        data_fname : str, required. 
            filename of input .fits data frame
        sky_fname : str, required. 
            filename of input .fits sky frame
        """

        super().__init__(data_fname, instrument)

        self.sky = Image(sky_fname).data

        # handle subarrays
        if self.sky.shape[0] != self.subc:
            self.sky = crop_center(self.sky, self.subc)

    def apply_sky(self):
        """
        Description
        -----------
        simply subtract the sky frame from the data
        """
        self.frames = np.array([data - self.sky for data in self.frames])

    def uranus_crop(self, bw):
        """
        Description
        -----------
        Custom crop for Uranus data. We use the right side of the
        NIRC2 detector to avoid the bad pixels in the lower left corner
        so just cutting off the left bit here
        
        Parameters
        ----------
        bw : int, required.
            width of the crop
        """
        szx, szy = self.frames[0].shape[0], self.frames[0].shape[1]
        self.frames = np.array([data[bw : szx - bw, 2 * bw :] for data in self.frames])


class Bxy3(Observation):
    """
    Description
    -----------
    The famous bxy3 dither at Keck,
    which avoids the noisier lower left quadrant of the detector
    input fnames MUST be in the conventional order 
    where target is in positions
    [top left, bottom right, top right]
        
    Workflow
    --------
    obs = bxy3.Bxy3(fnames)
    ## to check any step use:
    # obs.plot_frames()
    obs.make_sky(outdir+'sky_'+filt_name+'.fits')
    obs.apply_sky(outdir+'sky_'+filt_name+'.fits')
    obs.apply_flat(outdir+'flat_master_'+flat_filt+'.fits')
    obs.apply_badpx_map(outdir+'badpx_map_'+flat_filt+'.fits')
    obs.dewarp()
    obs.remove_cosmic_rays() # at this step there remain a few spots where pixel value is much lower than Neptune pixels. Doubt this is the fault of cosmic ray program; possibly failing to find them in flats with bad pixel search. Perhaps multi-layer search, e.g. large blocksize first, small blocksize second
    obs.rotate()
    obs.per_second()
    obs.write_frames([outdir+'frame0_nophot_'+filt_name+'.fits',outdir+'frame1_nophot_'+filt_name+'.fits',outdir+'frame2_nophot_'+filt_name+'.fits'])
    obs.trim()
    obs.stack()
    obs.crop_final(50)
    #check final
    plt.imshow(obs.final, origin = 'lower left')
    plt.show()
    obs.write_final(outdir+'stacked_nophot_'+filt_name+'.fits')#,png=True, png_file = outdir+target_name+'_'+filt_name+'.png')
    """

    def __init__(self, fnames, instrument):
        '''
        Parameters
        ----------
        fnames : list, required. 
            fits files representing a Keck bxy3 dither
        instrument : str, required.
        '''

        super().__init__(fnames, instrument)

    def make_sky(self, outfile):
        """
        Description
        -----------
        Make a sky frame via simple median-average of the frames
        Works because the planet is in only one quadrant at a time
        Writes fits file containing the sky frame
        with header info same as zeroth input bxy3 file
        
        Parameters
        ----------
        outfile : str, required.
            fname to write.
        """
        # make the master sky
        sky = np.median(self.frames, axis=0)
        # change some header info and write to .fits
        hdulist_out = self.dummy_fits.hdulist
        hdulist_out[0].header[self.header_kw_dict['object']] = "SKY_FULL"
        hdulist_out[0].data = sky
        hdulist_out[0].writeto(outfile, overwrite=True)

    def apply_sky(self, fname):
        """identify patches of "normal" sky in each of the three bxy3 frames. 
        comes up with a median average value of background for each frame. 
        Then you normalize the "full" sky frame we got from median of bxy3 
        to the value of the sky brightness in the single frame.
        Also normalizes everything to be around 1, so can directly multiply by data
        This ONLY works if bxy3 frames are input in the correct (usual) order,
        top left | bottom right | top right
        
        Parameters
        ----------
        fname : str, required. 
            filename pointing to sky frame
        """

        full_sky = Image(
            fname
        ).data  # could just as easily use self.sky, but this way could theoretically load a different sky background if desired
        if full_sky.shape[0] != self.subc:
            full_sky = crop_center(full_sky, self.subc)

        sz = int(self.subc / 16)  # just hardcoded this with a reasonable value
        c0 = (int(self.subc / 2 + self.subc / 4), int(self.subc / 2 - self.subc / 4))
        c1 = (int(self.subc / 2 - self.subc / 4), int(self.subc / 2 + self.subc / 4))
        c2 = (int(self.subc / 2 + self.subc / 4), int(self.subc / 2 + self.subc / 4))
        box0 = [c0[0] - sz, c0[0] + sz, c0[1] - sz, c0[1] + sz]
        box1 = [c1[0] - sz, c1[0] + sz, c1[1] - sz, c1[1] + sz]
        box2 = [c2[0] - sz, c2[0] + sz, c2[1] - sz, c2[1] + sz]
        # take bits of sky where planet is not
        skybits0 = [
            self.frames[0][box1[0] : box1[1], box1[2] : box1[3]],
            self.frames[0][box2[0] : box2[1], box2[2] : box2[3]],
        ]
        skybits1 = [
            self.frames[1][box0[0] : box0[1], box0[2] : box0[3]],
            self.frames[1][box2[0] : box2[1], box2[2] : box2[3]],
        ]
        skybits2 = [
            self.frames[2][box0[0] : box0[1], box0[2] : box0[3]],
            self.frames[2][box1[0] : box1[1], box1[2] : box1[3]],
        ]

        frames_skysub = []
        fullmed = np.median(full_sky)
        if np.abs(fullmed) > 5.0:  # minimum average counts hardcoded here
            for i in range(3):
                skybits = [skybits0, skybits1, skybits2][i]
                frame = self.frames[i]
                med0 = np.median(skybits[0])
                med1 = np.median(skybits[1])
                med = np.mean([med0, med1])
                norm = med / fullmed
                sky_norm = full_sky * norm
                frames_skysub.append(frame - sky_norm)
            self.frames = np.asarray(frames_skysub)
        else:
            print("Sky subtraction not applied (counts too low for good stats)")

    def trim(self):
        """
        Description
        -----------
        Clips each image in bxy3 to its own quadrant.
        Relies on frames input being in correct order.
        Should be applied just before stacking
        """
        sz = int(self.subc / 4)
        c0 = (int(self.subc / 2 + self.subc / 4), int(self.subc / 2 - self.subc / 4))
        c1 = (int(self.subc / 2 - self.subc / 4), int(self.subc / 2 + self.subc / 4))
        c2 = (int(self.subc / 2 + self.subc / 4), int(self.subc / 2 + self.subc / 4))
        box0 = [c0[0] - sz, c0[0] + sz, c0[1] - sz, c0[1] + sz]
        box1 = [c1[0] - sz, c1[0] + sz, c1[1] - sz, c1[1] + sz]
        box2 = [c2[0] - sz, c2[0] + sz, c2[1] - sz, c2[1] + sz]
        ftrim0 = self.frames[0][box0[0] : box0[1], box0[2] : box0[3]]
        ftrim1 = self.frames[1][box1[0] : box1[1], box1[2] : box1[3]]
        ftrim2 = self.frames[2][box2[0] : box2[1], box2[2] : box2[3]]
        self.frames = np.asarray([ftrim0, ftrim1, ftrim2])


class DitherN(Observation):
    '''
    Generic N-position dither; fnames can have arbitrary length
    will make sky by simple median average
    
    Workflow
    --------
    obs = DitherN(fnames, 'osiris')
    obs.make_sky(outdir+'sky_'+filt_name+'.fits')
    obs.apply_sky(outdir+'sky_'+filt_name+'.fits')
    obs.apply_flat(outdir+'flat_master_'+flat_filt+'.fits')
    obs.apply_badpx_map(outdir+'badpx_map_'+flat_filt+'.fits')
    obs.dewarp() #no distortion solution for osiris at time of writing
    obs.remove_cosmic_rays()
    obs.rotate() #no distortion solution for osiris at time of writing
    obs.per_second()
    obs.write_frames([outdir+'frame0_nophot_'+filt_name+'.fits',outdir+'frame1_nophot_'+filt_name+'.fits',outdir+'frame2_nophot_'+filt_name+'.fits'])
    obs.stack()
    obs.crop_final(500)
    obs.plot_final(show=True, png_file=outdir+'stacked_nophot_'+filt_name+'.png')
    obs.write_final(outdir+'stacked_nophot_'+filt_name+'.fits')
    '''
    
    def __init__(self, fnames, instrument):
        '''
        Parameters
        ----------
        fnames : list, required. 
            fits files representing a Keck bxy3 dither
        instrument : str, required.
        '''

        super().__init__(fnames, instrument)
        
    def make_sky(self, outfile):
        """
        Description
        -----------
        Make a sky frame via simple median-average of the frames
        Identical to Bxy3.make_sky()
        
        Parameters
        ----------
        outfile : str, required.
            fname to write.
        """
        # make the master sky
        sky = np.median(self.frames, axis=0)
        # change some header info and write to .fits
        hdulist_out = self.dummy_fits.hdulist
        hdulist_out[0].header[self.header_kw_dict['object']] = "SKY_FULL"
        hdulist_out[0].data = sky
        hdulist_out[0].writeto(outfile, overwrite=True)
    
    
    def apply_sky(self, fname):
        """
        Description
        -----------
        simply subtract the sky frame from the data
        
        Parameters
        ----------
        fname : str, required.
            filename representing sky frame fits
        """
        sky = Image(fname).data
        frames = np.array([data - sky for data in self.frames])
    