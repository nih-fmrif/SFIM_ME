import sys
import numpy as np
import nibabel as nib

__version__ = "0.10"

help_desc = """
This is a library that includes functions called for multi-echo (ME) automatic
denoising analyses. The functions included here are general utility functions
that aren't necessary specific to SFIM_ME. For example functions that help with
file loading or saving or module dependency checking
"""


def dep_check(ModuleSet):
    """
    This function takes a set of module names and makes sure all the
    modules are available to be imported.

    Parameters:
    -----------
    ModuleSet: Is a set of module names like set(["numpy", "pandas"])

    Results:
    --------
    Modules that aren't available are listed. If a module isn't available
        the program exits
    """

    try:
        import importlib
    except:
        import imp

    print("++ INFO [Main]: Checking for dependencies....")
    fails = 0

    for m in ModuleSet:
        try:
            modulefound = importlib.util.find_spec(m)
            if modulefound is None:
                print("++ ERROR [Main]: Can't import Module %s." +
                      " Please install." % m)
                fails += 1
        except:
            try:
                imp.find_module(m)
            except:
                print("++ ERROR [Main]: Can't import Module %s." +
                      " Please install." % m)
                fails += 1

    if fails == 0:
        print(" +              All Dependencies are OK.")
    else:
        print(" +               All dependencies not available. " +
              "Please install according to the above error messages.")
        print(" +               Program exited.")
        sys.exit()


def niiLoad(path):
    """
    This is a wrapper for the nibabel load function that loads the data
    and makes separate variables for the data, affine transform, and header
    """

    mepi_dset = nib.load(path)
    data = mepi_dset.get_data()
    aff = mepi_dset.affine
    head = mepi_dset.header
    # Note sure why it's necessary to change the
    # sform coed and null extensions here --DH
    head.extensions = []
    head.set_sform(head.get_sform(), code=1)
    return data, aff, head


def loadMEdata(data_file, Ne):
    """
    This function will load the multi-echo EPI datasets

    Parameters:
    -----------
    data_file: this can be a single file name or a comma delimited list of
      files. If it's a single file name, the file is a 4D data set where
      the multiple echos are concatinated in the z direction with the
      shortest echo having the lowest z indices.
      If this is a list of file names, each file is a 4D data set from
      a single echo starting with the shortest echo.
    Ne: Number of echos

    Returns:
    --------
    mepi_dset_data: A matrix with 5D data that's organized as:
                    X,Y,Z,NumEchos,NumTimePoints
    Nx, Ny, Nz, Nt: The # of values in the X, Y, Z, and time dimensions
    mepi_deset_aff: The affine transform for the data set
    mepi_dset_head: The header from the dataset with the extensions
                    subfield nulled

    If multiple files are inputted and Nx, Ny, Nz, Nt, and mepi_deset_aff
    don't match across volumes, this function will give an error and not
    output datasets
    """

    # Loading data from one file with the echos concatinated in the z direction
    # -----------------------------------------------------------------------
    if ',' not in data_file:
        print("++ INFO: Attemting to load input dataset %s" % data_file)
        mepi_dset_data, mepi_dset_aff, mepi_dset_head = niiLoad(data_file)
        Ndim = mepi_dset_data.ndim
        if Ndim != 4:
            print("++ Error: The input dataset [%s] has the wrong number of " +
                  "dimensions [Ndim=%i]." % (data_file, Ndim))
            sys.exit()

        # Control correct size of dataset on the Z-direction
        # --------------------------------------------------
        Nx, Ny, Nz, Nt = mepi_dset_data.shape
        if (Nz % Ne == 0):
            Nz = int(Nz/Ne)
        else:
            print("++ Error: The Z-dimension in dataset [%s] is not a " +
                  "multiple of the number of echo times [Ne=%i]."
                  % (data_file, Ne))
            print("++ INFO: Input Dataset Dimensions: " +
                  "[Nx=%i,Ne=%i,Nz=%i,Ne=%i,Nt=%i]." % (Nx, Ny, Nz, Ne, Nt))
            sys.exit()

        # Reshaping the 4D data set with concatinated echos to a 5D data set
        mepi_dset_data = mepi_dset_data.reshape((Nx, Ny, Nz, Ne, Nt),
                                                order='F')
    # Loading data from a comma separated list of files with each being a
    # separate echo
    # -----------------------------------------------------------------------
    else:
        data_files = data_file.split(',')
        print("--------------------------------------------------")
        print("++INFO: Attemting to load input from multiple datasets:")
        for i in range(len(data_files)):
            print(data_files[i])
        print("--------------------------------------------------")
        if (Ne != len(data_files)):
            print("++Error: There are " + Ne + " echo times but only " +
                  len(data_files) + " file names")
            sys.exit()

        for i in range(Ne):
            mepi_dset_data_tmp, mepi_dset_aff_tmp, mepi_dset_head_tmp = \
                niiLoad(data_files[i])
            Ndim_tmp = mepi_dset_data_tmp.ndim
            if Ndim_tmp != 4:
                print("++ Error: The input dataset [%s] has the wrong number" +
                      " of dimensions [Ndim=%i]." % (data_files[i], Ndim))
                sys.exit()

            if i == 0:
                # Initalizing volume parameters for the first echo volume
                Nx, Ny, Nz, Nt = mepi_dset_data_tmp.shape
                mepi_dset_data = np.empty((Nx, Ny, Nz, Ne, Nt))
                mepi_dset_data[:, :, :, i, :] = mepi_dset_data_tmp
                mepi_dset_aff = mepi_dset_aff_tmp
                mepi_dset_head = mepi_dset_head_tmp
                pixdim = mepi_dset_head['pixdim']
            else:
                # For echo volumes after the first one is initalized
                # Check if the affine matrix is the same between the first
                # and current echo
                affcomptmp = mepi_dset_aff == mepi_dset_aff_tmp
                if not affcomptmp.all():
                    print("++ Error: "+data_files[0] + " and " + data_files[1]
                          + " have different affine transform matrices")
                    print(mepi_dset_aff)
                    print("versus")
                    print(mepi_dset_aff_tmp)
                    sys.exit()

                # Check if the voxel grid is the same between the first and
                # current echo
                Nxtmp, Nytmp, Nztmp, Nttmp = mepi_dset_data_tmp.shape
                if not ((Nx == Nxtmp) & (Ny == Nytmp) &
                        (Nz == Nztmp) & (Nt == Nttmp)):
                    print("++ Error: "+data_files[0] + " and " + data_files[1]
                          + " have different matrix sizes")
                    print("(X,Y,Z,T) is (%d, %d, %d, %d) vs. (%d, %d, %d, %d)"
                          % (Nx, Ny, Nz, Nt, Nxtmp, Nytmp, Nztmp, Nttmp))
                    sys.exit()

                # Check to make sure voxel sizes are the same across echos
                pixdimtmp = mepi_dset_head_tmp['pixdim']
                pixdimcomptmp = pixdim[1:3] == pixdimtmp[1:3]
                if not pixdimcomptmp.all():
                    print("++ Error: "+data_files[0] + " and " + data_files[1]
                          + " have different voxel sizes")
                    print(pixdim[1:3] + "versus" + pixdimtmp[1:3])
                    sys.exit()

                # ANYTHING ELSE TO TEST HERE???
                # Make sure TRs match in all files???

                # Assuming everything else is the same, the only parameter
                # that needs to be added per volume is the voxel data
                mepi_dset_data[:, :, :, i] = mepi_dset_data_tmp

    return mepi_dset_data, Nx, Ny, Nz, Ny, Nt, mepi_dset_aff, mepi_dset_head


def niiwrite_nv(data, mask, temp_path, aff, temp_header):
    """
    This function will write NIFTI datasets

    Parameters:
    ----------
    data: this is (Nv, Nt) or (Nv,) array. No z-cat ME datasets allowed.
    mask: this is (Nx,Ny,Nz) array with Nv entries equal to True.
          This is an intracranial voxel mask.
    temp_path: this is the output directory.
    aff: affine transformation associated with this dataset.
    temp_header: header for the dataset.

    Returns:
    --------
    None.
    """
    Nx, Ny, Nz = mask.shape
    if (data.ndim == 1):
        temp = np.zeros((Nx, Ny, Nz), order='F')
        temp[mask] = data
    if (data.ndim == 2):
            _, Nt = data.shape
            temp = np.zeros((Nx, Ny, Nz, Nt), order='F')
            temp[mask, :] = data
    if (data.ndim == 3):
            Nv, Ne, Nt = data.shape
            temp = np.zeros((Nx, Ny, Nz, Nt), order='F')
            temp[mask, :] = data[:, 0, :]
            for e in range(1, Ne):
                aux = np.zeros((Nx, Ny, Nz, Nt), order='F')
                aux[mask, :] = data[:, e, :]
                temp = np.concatenate((temp, aux), axis=2)

    outni = nib.Nifti1Image(temp, aff, header=temp_header)
    outni.to_filename(temp_path)
    print(" +              Dataset %s written to disk" % (temp_path))


def mask4MEdata(data):
    """
    This function will create a mask for ME data taking into
      account the time series of all echoes

    Paramters:
    ----------
    data: this is a (Nx,Ny,Nz,Ne,Nt) array with ME data.

    Returns:
    --------
    mask: this is a (Nx,Ny,Nz) marking all voxels that have valid time-series
          for all echo times.
    """
    # Create mask taking into account all echoes
    Nx, Ny, Nz, Ne, Nt = data.shape
    mask = np.ones((Nx, Ny, Nz), dtype=np.bool)
    for i in range(Ne):
        tmpmask = (data[:, :, :, i, :] != 0).prod(axis=-1, dtype=np.bool)
        mask = mask & tmpmask
    return mask
