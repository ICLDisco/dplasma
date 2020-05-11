/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 */
#ifndef _DPLASMA_LU_QR_H_
#define _DPLASMA_LU_QR_H_

BEGIN_C_DECLS

/*
 * Enum criteria for LU/QR algorithm
 */
enum dplasma_criteria_e {
    DEFAULT_CRITERIUM    = 0,
    HIGHAM_CRITERIUM     = 1,
    MUMPS_CRITERIUM      = 2,
    LU_ONLY_CRITERIUM    = 3,
    QR_ONLY_CRITERIUM    = 4,
    RANDOM_CRITERIUM     = 5,
    HIGHAM_SUM_CRITERIUM = 6,
    HIGHAM_MAX_CRITERIUM = 7,
    HIGHAM_MOY_CRITERIUM = 8
};

END_C_DECLS

#endif /* _DPLASMA_LU_QR_H_ */
