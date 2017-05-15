

MAT_CHOICES = {1: 'accept', 2: 'maybe', 3: 'reject'}


def load_mat_validation(mat_paths):
    # scipy breaks Ctrl-C interrupts when scipy for some distributions,
    # only import as necessary
    from scipy.io import loadmat

    validation_data = {}

    for mat_path in mat_paths:
        mat_data = loadmat(mat_path)

        for scan in mat_data['data'][0]:
            scan_num = scan['scan_number'][0][0][0][0]

            if 'fragments' not in scan.dtype.names:
                continue

            for ptm_combo in scan['fragments'][0][0][0]:
                seq = ptm_combo['seq'][0][0][0]
                choice = ptm_combo['status'][0][0][0][0]

                if choice == 0:
                    continue

                validation_data[scan_num, seq] = MAT_CHOICES[choice]

    return validation_data
