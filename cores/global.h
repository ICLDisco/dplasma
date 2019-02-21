/**
 *
 * @file global.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @date 2010-11-15
 *
 **/

/***************************************************************************//**
 *  PLASMA internals of interest to PLASMA core developers, but not necessarily
 *  of interest to PLASMA community contributors.
 **/
#ifndef _PLASMA_GLOBAL_H_
#define _PLASMA_GLOBAL_H_

#include <dplasma.h>

#include <assert.h>
#include <string.h>

/***************************************************************************//**
 *  Configuration
 **/
// maximum contexts
#define CONTEXTS_MAX         256
// maximum cores per context
#define CONTEXT_THREADS_MAX  256
// size of parallel functions arguments buffer
#define ARGS_BUFF_SIZE       512
// cache line size
#define CACHE_LINE_SIZE      128
// standard page size
#define STANDARD_PAGE_SIZE  4096

/***************************************************************************//**
 *  Action commands
 **/
#define PLASMA_ACT_STAND_BY     0
#define PLASMA_ACT_PARALLEL     1
#define PLASMA_ACT_DYNAMIC      2
#define PLASMA_ACT_FINALIZE     3

/***************************************************************************//**
 *  Numerical operations
 **/
#define PLASMA_FUNC_SGELS    1
#define PLASMA_FUNC_SPOSV    2
#define PLASMA_FUNC_SGESV    3
#define PLASMA_FUNC_DGELS    4
#define PLASMA_FUNC_DPOSV    5
#define PLASMA_FUNC_DGESV    6
#define PLASMA_FUNC_CGELS    7
#define PLASMA_FUNC_CPOSV    8
#define PLASMA_FUNC_CGESV    9
#define PLASMA_FUNC_ZGELS   10
#define PLASMA_FUNC_ZPOSV   11
#define PLASMA_FUNC_ZGESV   12
#define PLASMA_FUNC_ZCGESV  13
#define PLASMA_FUNC_DSGESV  14
#define PLASMA_FUNC_ZCPOSV  15
#define PLASMA_FUNC_DSPOSV  16
#define PLASMA_FUNC_DSGELS  17
#define PLASMA_FUNC_ZCGELS  18
#define PLASMA_FUNC_SGEMM   19
#define PLASMA_FUNC_DGEMM   20
#define PLASMA_FUNC_CGEMM   21
#define PLASMA_FUNC_ZGEMM   22
#define PLASMA_FUNC_SSYMM   23
#define PLASMA_FUNC_DSYMM   24
#define PLASMA_FUNC_CSYMM   25
#define PLASMA_FUNC_ZSYMM   26
#define PLASMA_FUNC_CHERK   27
#define PLASMA_FUNC_ZHERK   28
#define PLASMA_FUNC_SSYRK   29
#define PLASMA_FUNC_DSYRK   30
#define PLASMA_FUNC_CSYRK   31
#define PLASMA_FUNC_ZSYRK   32
#define PLASMA_FUNC_CHEMM   33
#define PLASMA_FUNC_ZHEMM   34
#define PLASMA_FUNC_ZHEEV   35
#define PLASMA_FUNC_CHEEV   36
#define PLASMA_FUNC_DSYEV   37
#define PLASMA_FUNC_SSYEV   38
#define PLASMA_FUNC_ZHEEVD  39
#define PLASMA_FUNC_CHEEVD  40
#define PLASMA_FUNC_DSYEVD  41
#define PLASMA_FUNC_SSYEVD  42
#define PLASMA_FUNC_ZHEGST  43
#define PLASMA_FUNC_CHEGST  44
#define PLASMA_FUNC_DSYGST  45
#define PLASMA_FUNC_SSYGST  46
#define PLASMA_FUNC_ZHEGV   47
#define PLASMA_FUNC_CHEGV   48
#define PLASMA_FUNC_DSYGV   49
#define PLASMA_FUNC_SSYGV   50
#define PLASMA_FUNC_ZHEGVD  51
#define PLASMA_FUNC_CHEGVD  52
#define PLASMA_FUNC_DSYGVD  53
#define PLASMA_FUNC_SSYGVD  54
#define PLASMA_FUNC_ZHETRD  55
#define PLASMA_FUNC_CHETRD  56
#define PLASMA_FUNC_DSYTRD  57
#define PLASMA_FUNC_SSYTRD  58
#define PLASMA_FUNC_ZGESVD  59
#define PLASMA_FUNC_CGESVD  60
#define PLASMA_FUNC_DGESVD  61
#define PLASMA_FUNC_SGESVD  62
#define PLASMA_FUNC_ZGEEV   63
#define PLASMA_FUNC_CGEEV   64
#define PLASMA_FUNC_DGEEV   65
#define PLASMA_FUNC_SGEEV   66
#define PLASMA_FUNC_ZGEHRD  67
#define PLASMA_FUNC_CGEHRD  68
#define PLASMA_FUNC_DGEHRD  69
#define PLASMA_FUNC_SGEHRD  70
#define PLASMA_FUNC_ZGEBRD  71
#define PLASMA_FUNC_CGEBRD  72
#define PLASMA_FUNC_DGEBRD  73
#define PLASMA_FUNC_SGEBRD  74

/***************************************************************************//**
 *  Parallel function call - packing of arguments
 **/
#define plasma_pack_args_1( \
    type1, arg1) \
{ \
    type1 var1 = (arg1); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_1", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
}

#define plasma_pack_args_2( \
    type1, arg1, \
    type2, arg2) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_2", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
}

#define plasma_pack_args_3( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_3", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
}

#define plasma_pack_args_4( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_4", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
}

#define plasma_pack_args_5( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_5", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
}

#define plasma_pack_args_6( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_6", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
}

#define plasma_pack_args_7( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_7", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
}

#define plasma_pack_args_8( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_8", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
}

#define plasma_pack_args_9( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_9", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
}

#define plasma_pack_args_10( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    type10 var10 = (arg10); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) + \
        sizeof(type10) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_9", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
    memcpy(plasma_ptr, &var10, sizeof(type10)); plasma_ptr += sizeof(type10); \
}


#define plasma_pack_args_11( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    type10 var10 = (arg10); \
    type11 var11 = (arg11); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) + \
        sizeof(type10) + \
        sizeof(type11) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_11", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
    memcpy(plasma_ptr, &var10, sizeof(type10)); plasma_ptr += sizeof(type10); \
    memcpy(plasma_ptr, &var11, sizeof(type11)); plasma_ptr += sizeof(type11); \
}

#define plasma_pack_args_12( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    type10 var10 = (arg10); \
    type11 var11 = (arg11); \
    type12 var12 = (arg12); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) + \
        sizeof(type10) + \
        sizeof(type11) + \
        sizeof(type12) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_12", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
    memcpy(plasma_ptr, &var10, sizeof(type10)); plasma_ptr += sizeof(type10); \
    memcpy(plasma_ptr, &var11, sizeof(type11)); plasma_ptr += sizeof(type11); \
    memcpy(plasma_ptr, &var12, sizeof(type12)); plasma_ptr += sizeof(type12); \
}


#define plasma_pack_args_13( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    type10 var10 = (arg10); \
    type11 var11 = (arg11); \
    type12 var12 = (arg12); \
    type13 var13 = (arg13); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) + \
        sizeof(type10) + \
        sizeof(type11) + \
        sizeof(type12) + \
        sizeof(type13) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_13", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
    memcpy(plasma_ptr, &var10, sizeof(type10)); plasma_ptr += sizeof(type10); \
    memcpy(plasma_ptr, &var11, sizeof(type11)); plasma_ptr += sizeof(type11); \
    memcpy(plasma_ptr, &var12, sizeof(type12)); plasma_ptr += sizeof(type12); \
    memcpy(plasma_ptr, &var13, sizeof(type13)); plasma_ptr += sizeof(type13); \
}



#define plasma_pack_args_14( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    type10 var10 = (arg10); \
    type11 var11 = (arg11); \
    type12 var12 = (arg12); \
    type13 var13 = (arg13); \
    type14 var14 = (arg14); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) + \
        sizeof(type10) + \
        sizeof(type11) + \
        sizeof(type12) + \
        sizeof(type13) + \
        sizeof(type14) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_14", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
    memcpy(plasma_ptr, &var10, sizeof(type10)); plasma_ptr += sizeof(type10); \
    memcpy(plasma_ptr, &var11, sizeof(type11)); plasma_ptr += sizeof(type11); \
    memcpy(plasma_ptr, &var12, sizeof(type12)); plasma_ptr += sizeof(type12); \
    memcpy(plasma_ptr, &var13, sizeof(type13)); plasma_ptr += sizeof(type13); \
    memcpy(plasma_ptr, &var14, sizeof(type14)); plasma_ptr += sizeof(type14); \
}

#define plasma_pack_args_15( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14, \
    type15, arg15) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    type10 var10 = (arg10); \
    type11 var11 = (arg11); \
    type12 var12 = (arg12); \
    type13 var13 = (arg13); \
    type14 var14 = (arg14); \
    type15 var15 = (arg15); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) + \
        sizeof(type10) + \
        sizeof(type11) + \
        sizeof(type12) + \
        sizeof(type13) + \
        sizeof(type14) + \
        sizeof(type15) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_15", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
    memcpy(plasma_ptr, &var10, sizeof(type10)); plasma_ptr += sizeof(type10); \
    memcpy(plasma_ptr, &var11, sizeof(type11)); plasma_ptr += sizeof(type11); \
    memcpy(plasma_ptr, &var12, sizeof(type12)); plasma_ptr += sizeof(type12); \
    memcpy(plasma_ptr, &var13, sizeof(type13)); plasma_ptr += sizeof(type13); \
    memcpy(plasma_ptr, &var14, sizeof(type14)); plasma_ptr += sizeof(type14); \
    memcpy(plasma_ptr, &var15, sizeof(type15)); plasma_ptr += sizeof(type15); \
}


#define plasma_pack_args_16( \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14, \
    type15, arg15, \
    type16, arg16) \
{ \
    type1 var1 = (arg1); \
    type2 var2 = (arg2); \
    type3 var3 = (arg3); \
    type4 var4 = (arg4); \
    type5 var5 = (arg5); \
    type6 var6 = (arg6); \
    type7 var7 = (arg7); \
    type8 var8 = (arg8); \
    type9 var9 = (arg9); \
    type10 var10 = (arg10); \
    type11 var11 = (arg11); \
    type12 var12 = (arg12); \
    type13 var13 = (arg13); \
    type14 var14 = (arg14); \
    type15 var15 = (arg15); \
    type16 var16 = (arg16); \
    unsigned char *plasma_ptr = plasma->args_buff; \
    if (sizeof(type1) + \
        sizeof(type2) + \
        sizeof(type3) + \
        sizeof(type4) + \
        sizeof(type5) + \
        sizeof(type6) + \
        sizeof(type7) + \
        sizeof(type8) + \
        sizeof(type9) + \
        sizeof(type10) + \
        sizeof(type11) + \
        sizeof(type12) + \
        sizeof(type13) + \
        sizeof(type14) + \
        sizeof(type15) + \
        sizeof(type16) > ARGS_BUFF_SIZE) \
        plasma_fatal_error("plasma_pack_args_16", "arguments buffer too small"); \
    memcpy(plasma_ptr, &var1, sizeof(type1)); plasma_ptr += sizeof(type1); \
    memcpy(plasma_ptr, &var2, sizeof(type2)); plasma_ptr += sizeof(type2); \
    memcpy(plasma_ptr, &var3, sizeof(type3)); plasma_ptr += sizeof(type3); \
    memcpy(plasma_ptr, &var4, sizeof(type4)); plasma_ptr += sizeof(type4); \
    memcpy(plasma_ptr, &var5, sizeof(type5)); plasma_ptr += sizeof(type5); \
    memcpy(plasma_ptr, &var6, sizeof(type6)); plasma_ptr += sizeof(type6); \
    memcpy(plasma_ptr, &var7, sizeof(type7)); plasma_ptr += sizeof(type7); \
    memcpy(plasma_ptr, &var8, sizeof(type8)); plasma_ptr += sizeof(type8); \
    memcpy(plasma_ptr, &var9, sizeof(type9)); plasma_ptr += sizeof(type9); \
    memcpy(plasma_ptr, &var10, sizeof(type10)); plasma_ptr += sizeof(type10); \
    memcpy(plasma_ptr, &var11, sizeof(type11)); plasma_ptr += sizeof(type11); \
    memcpy(plasma_ptr, &var12, sizeof(type12)); plasma_ptr += sizeof(type12); \
    memcpy(plasma_ptr, &var13, sizeof(type13)); plasma_ptr += sizeof(type13); \
    memcpy(plasma_ptr, &var14, sizeof(type14)); plasma_ptr += sizeof(type14); \
    memcpy(plasma_ptr, &var15, sizeof(type15)); plasma_ptr += sizeof(type15); \
    memcpy(plasma_ptr, &var16, sizeof(type16)); plasma_ptr += sizeof(type16); \
}

/***************************************************************************//**
 *  Sync after dynamically scheduled section
 **/
#define plasma_dynamic_sync() \
{ \
    if (plasma->dynamic_section) { \
        QUARK_Waitall(plasma->quark); \
        plasma_barrier(plasma); \
        plasma->dynamic_section = PLASMA_FALSE; \
    } \
}

/***************************************************************************//**
 *  Parallel SPMD function call - thread control
 **/
#define plasma_static_call(parallel_function) \
{ \
    if (plasma->dynamic_section) \
        plasma_dynamic_sync(); \
    pthread_mutex_lock(&plasma->action_mutex); \
    plasma->action = PLASMA_ACT_PARALLEL; \
    plasma->parallel_func_ptr = &parallel_function; \
    pthread_mutex_unlock(&plasma->action_mutex); \
    pthread_cond_broadcast(&plasma->action_condt); \
    plasma_barrier(plasma); \
    plasma->action = PLASMA_ACT_STAND_BY; \
    parallel_function(plasma); \
    plasma_barrier(plasma); \
}

/***************************************************************************//**
 *  Start dynamically scheduled section
 **/
#define plasma_dynamic_spawn() \
{ \
    if (!plasma->dynamic_section) { \
        plasma->dynamic_section = PLASMA_TRUE; \
        pthread_mutex_lock(&plasma->action_mutex); \
        plasma->action = PLASMA_ACT_DYNAMIC; \
        pthread_mutex_unlock(&plasma->action_mutex); \
        pthread_cond_broadcast(&plasma->action_condt); \
        plasma_barrier(plasma); \
        plasma->action = PLASMA_ACT_STAND_BY; \
    } \
}

/***************************************************************************//**
 *  Parallel call for functions with static versions only
 **/
#define plasma_static_call_1( \
           parallel_function, \
    type1, arg1) \
    plasma_pack_args_1( \
        type1, (arg1)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_2( \
           parallel_function, \
    type1, arg1, \
    type2, arg2) \
    plasma_pack_args_2( \
        type1, (arg1), \
        type2, (arg2)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_3( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3) \
    plasma_pack_args_3( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_4( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4) \
    plasma_pack_args_4( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_5( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5) \
    plasma_pack_args_5( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_6( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6) \
    plasma_pack_args_6( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_7( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7) \
    plasma_pack_args_7( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_8( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8) \
    plasma_pack_args_8( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_9( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9) \
    plasma_pack_args_9( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8), \
        type9, (arg9)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_10( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10) \
    plasma_pack_args_10( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8), \
        type9, (arg9), \
        type10, (arg10)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_11( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11) \
    plasma_pack_args_11( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8), \
        type9, (arg9), \
        type10, (arg10), \
        type11, (arg11)) \
    plasma_static_call(parallel_function) \


#define plasma_static_call_12( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12) \
    plasma_pack_args_12( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8), \
        type9, (arg9), \
        type10, (arg10), \
        type11, (arg11), \
        type12, (arg12)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_13( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13) \
    plasma_pack_args_13( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8), \
        type9, (arg9), \
        type10, (arg10), \
        type11, (arg11), \
        type12, (arg12), \
        type13, (arg13)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_14( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14) \
    plasma_pack_args_14( \
        type1, (arg1), \
        type2, (arg2), \
        type3, (arg3), \
        type4, (arg4), \
        type5, (arg5), \
        type6, (arg6), \
        type7, (arg7), \
        type8, (arg8), \
        type9, (arg9), \
        type10, (arg10), \
        type11, (arg11), \
        type12, (arg12), \
        type13, (arg13), \
        type14, (arg14)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_15( \
                           parallel_function, \
                    type1, arg1, \
                    type2, arg2, \
                    type3, arg3, \
                    type4, arg4, \
                    type5, arg5, \
                    type6, arg6, \
                    type7, arg7, \
                    type8, arg8, \
                    type9, arg9, \
                    type10, arg10, \
                    type11, arg11, \
                    type12, arg12, \
                    type13, arg13, \
                    type14, arg14, \
                    type15, arg15) \
    plasma_pack_args_15( \
                            type1, (arg1), \
                            type2, (arg2), \
                            type3, (arg3), \
                            type4, (arg4), \
                            type5, (arg5), \
                            type6, (arg6), \
                            type7, (arg7), \
                            type8, (arg8), \
                            type9, (arg9), \
                            type10, (arg10), \
                            type11, (arg11), \
                            type12, (arg12), \
                            type13, (arg13), \
                            type14, (arg14), \
                            type15, (arg15)) \
    plasma_static_call(parallel_function) \

#define plasma_static_call_16( \
                           parallel_function, \
                    type1, arg1, \
                    type2, arg2, \
                    type3, arg3, \
                    type4, arg4, \
                    type5, arg5, \
                    type6, arg6, \
                    type7, arg7, \
                    type8, arg8, \
                    type9, arg9, \
                    type10, arg10, \
                    type11, arg11, \
                    type12, arg12, \
                    type13, arg13, \
                    type14, arg14, \
                    type15, arg15, \
                    type16, arg16) \
    plasma_pack_args_16( \
                            type1, (arg1), \
                            type2, (arg2), \
                            type3, (arg3), \
                            type4, (arg4), \
                            type5, (arg5), \
                            type6, (arg6), \
                            type7, (arg7), \
                            type8, (arg8), \
                            type9, (arg9), \
                            type10, (arg10), \
                            type11, (arg11), \
                            type12, (arg12), \
                            type13, (arg13), \
                            type14, (arg14), \
                            type15, (arg15), \
                            type16, (arg16)) \
    plasma_static_call(parallel_function) \

/***************************************************************************//**
 *  Parallel call for functions with both static and dynamic versions
 **/
#define plasma_parallel_call_1( \
           parallel_function, \
    type1, arg1) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_1( \
            type1, (arg1)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1); \
    }

#define plasma_parallel_call_2( \
           parallel_function, \
    type1, arg1, \
    type2, arg2) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_2( \
            type1, (arg1), \
            type2, (arg2)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2); \
    }

#define plasma_parallel_call_3( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_3( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3); \
    }

#define plasma_parallel_call_4( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_4( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4); \
    }

#define plasma_parallel_call_5( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_5( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5); \
    }

#define plasma_parallel_call_6( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_6( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6); \
    }

#define plasma_parallel_call_7( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_7( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7); \
    }

#define plasma_parallel_call_8( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_8( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8); \
    }

#define plasma_parallel_call_9( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_9( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9); \
    }

#define plasma_parallel_call_10( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_10( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9), \
            type10, (arg10)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9, \
            arg10); \
    }

#define plasma_parallel_call_11( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_11( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9), \
            type10, (arg10), \
            type11, (arg11)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9, \
            arg10, \
            arg11); \
    }



#define plasma_parallel_call_12( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_12( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9), \
            type10, (arg10), \
            type11, (arg11), \
            type12, (arg12)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9, \
            arg10, \
            arg11, \
            arg12); \
    }


#define plasma_parallel_call_13( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_13( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9), \
            type10, (arg10), \
            type11, (arg11), \
            type12, (arg12), \
            type13, (arg13)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9, \
            arg10, \
            arg11, \
            arg12, \
            arg13); \
    }

#define plasma_parallel_call_14( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_14( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9), \
            type10, (arg10), \
            type11, (arg11), \
            type12, (arg12), \
            type13, (arg13), \
            type14, (arg14)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9, \
            arg10, \
            arg11, \
            arg12, \
            arg13, \
            arg14); \
    }
#define plasma_parallel_call_15( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14, \
    type15, arg15) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_15( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9), \
            type10, (arg10), \
            type11, (arg11), \
            type12, (arg12), \
            type13, (arg13), \
            type14, (arg14), \
            type15, (arg15)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9, \
            arg10, \
            arg11, \
            arg12, \
            arg13, \
            arg14, \
            arg15); \
    }
#define plasma_parallel_call_16( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14, \
    type15, arg15, \
    type16, arg16) \
    if (PLASMA_SCHEDULING == PLASMA_STATIC_SCHEDULING) { \
        plasma_pack_args_16( \
            type1, (arg1), \
            type2, (arg2), \
            type3, (arg3), \
            type4, (arg4), \
            type5, (arg5), \
            type6, (arg6), \
            type7, (arg7), \
            type8, (arg8), \
            type9, (arg9), \
            type10, (arg10), \
            type11, (arg11), \
            type12, (arg12), \
            type13, (arg13), \
            type14, (arg14), \
            type15, (arg15), \
            type16, (arg16)) \
        plasma_static_call(parallel_function) \
    } else { \
        plasma_dynamic_spawn(); \
        parallel_function##_quark( \
            arg1, \
            arg2, \
            arg3, \
            arg4, \
            arg5, \
            arg6, \
            arg7, \
            arg8, \
            arg9, \
            arg10, \
            arg11, \
            arg12, \
            arg13, \
            arg14, \
            arg15, \
            arg16); \
    }

/***************************************************************************//**
 *  Parallel call for functions with dynamic versions only
 **/
#define plasma_dynamic_call_1( \
           parallel_function, \
    type1, arg1) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1); \

#define plasma_dynamic_call_2( \
           parallel_function, \
    type1, arg1, \
    type2, arg2) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2); \

#define plasma_dynamic_call_3( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3); \

#define plasma_dynamic_call_4( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4); \

#define plasma_dynamic_call_5( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5); \

#define plasma_dynamic_call_6( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6); \

#define plasma_dynamic_call_7( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7); \

#define plasma_dynamic_call_8( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8); \

#define plasma_dynamic_call_9( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9);

#define plasma_dynamic_call_10( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9, \
        arg10);


#define plasma_dynamic_call_11( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9, \
        arg10, \
        arg11);


#define plasma_dynamic_call_12( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9, \
        arg10, \
        arg11, \
        arg12);

#define plasma_dynamic_call_13( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9, \
        arg10, \
        arg11, \
        arg12, \
        arg13);

#define plasma_dynamic_call_14( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9, \
        arg10, \
        arg11, \
        arg12, \
        arg13, \
        arg14);

#define plasma_dynamic_call_15( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14, \
    type15, arg15) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9, \
        arg10, \
        arg11, \
        arg12, \
        arg13, \
        arg14, \
        arg15);

#define plasma_dynamic_call_16( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14, \
    type15, arg15, \
    type16, arg16) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9, \
        arg10, \
        arg11, \
        arg12, \
        arg13, \
        arg14, \
        arg15, \
        arg16);

#define plasma_dynamic_call_21( \
           parallel_function, \
    type1, arg1, \
    type2, arg2, \
    type3, arg3, \
    type4, arg4, \
    type5, arg5, \
    type6, arg6, \
    type7, arg7, \
    type8, arg8, \
    type9, arg9, \
    type10, arg10, \
    type11, arg11, \
    type12, arg12, \
    type13, arg13, \
    type14, arg14, \
    type15, arg15, \
    type16, arg16, \
    type17, arg17, \
    type18, arg18, \
    type19, arg19, \
    type20, arg20, \
    type21, arg21) \
    plasma_dynamic_spawn(); \
    parallel_function##_quark( \
        arg1, \
        arg2, \
        arg3, \
        arg4, \
        arg5, \
        arg6, \
        arg7, \
        arg8, \
        arg9, \
        arg10, \
        arg11, \
        arg12, \
        arg13, \
        arg14, \
        arg15, \
        arg16, \
        arg17, \
        arg18, \
        arg19, \
        arg20, \
        arg21);
/***************************************************************************//**
 *  Parallel function call - unpacking of arguments
 **/
#define plasma_unpack_args_1( \
    arg1) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
}

#define plasma_unpack_args_2( \
    arg1, \
    arg2) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
}

#define plasma_unpack_args_3( \
    arg1, \
    arg2, \
    arg3) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
}

#define plasma_unpack_args_4( \
    arg1, \
    arg2, \
    arg3, \
    arg4) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
}

#define plasma_unpack_args_5( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
}

#define plasma_unpack_args_6( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
}

#define plasma_unpack_args_7( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
}

#define plasma_unpack_args_8( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
}

#define plasma_unpack_args_9( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
}

#define plasma_unpack_args_10( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
}

#define plasma_unpack_args_11( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10, \
    arg11) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
    memcpy(&arg11, plasma_ptr, sizeof(arg11)); plasma_ptr += sizeof(arg11); \
}


#define plasma_unpack_args_12( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10, \
    arg11, \
    arg12) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
    memcpy(&arg11, plasma_ptr, sizeof(arg11)); plasma_ptr += sizeof(arg11); \
    memcpy(&arg12, plasma_ptr, sizeof(arg12)); plasma_ptr += sizeof(arg12); \
}

#define plasma_unpack_args_13( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10, \
    arg11, \
    arg12, \
    arg13) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
    memcpy(&arg11, plasma_ptr, sizeof(arg11)); plasma_ptr += sizeof(arg11); \
    memcpy(&arg12, plasma_ptr, sizeof(arg12)); plasma_ptr += sizeof(arg12); \
    memcpy(&arg13, plasma_ptr, sizeof(arg13)); plasma_ptr += sizeof(arg13); \
}

#define plasma_unpack_args_14( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10, \
    arg11, \
    arg12, \
    arg13, \
    arg14) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
    memcpy(&arg11, plasma_ptr, sizeof(arg11)); plasma_ptr += sizeof(arg11); \
    memcpy(&arg12, plasma_ptr, sizeof(arg12)); plasma_ptr += sizeof(arg12); \
    memcpy(&arg13, plasma_ptr, sizeof(arg13)); plasma_ptr += sizeof(arg13); \
    memcpy(&arg14, plasma_ptr, sizeof(arg14)); plasma_ptr += sizeof(arg14); \
}

#define plasma_unpack_args_15( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10, \
    arg11, \
    arg12, \
    arg13, \
    arg14, \
    arg15) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
    memcpy(&arg11, plasma_ptr, sizeof(arg11)); plasma_ptr += sizeof(arg11); \
    memcpy(&arg12, plasma_ptr, sizeof(arg12)); plasma_ptr += sizeof(arg12); \
    memcpy(&arg13, plasma_ptr, sizeof(arg13)); plasma_ptr += sizeof(arg13); \
    memcpy(&arg14, plasma_ptr, sizeof(arg14)); plasma_ptr += sizeof(arg14); \
    memcpy(&arg15, plasma_ptr, sizeof(arg15)); plasma_ptr += sizeof(arg15); \
}

#define plasma_unpack_args_16( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10, \
    arg11, \
    arg12, \
    arg13, \
    arg14, \
    arg15, \
    arg16) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
    memcpy(&arg11, plasma_ptr, sizeof(arg11)); plasma_ptr += sizeof(arg11); \
    memcpy(&arg12, plasma_ptr, sizeof(arg12)); plasma_ptr += sizeof(arg12); \
    memcpy(&arg13, plasma_ptr, sizeof(arg13)); plasma_ptr += sizeof(arg13); \
    memcpy(&arg14, plasma_ptr, sizeof(arg14)); plasma_ptr += sizeof(arg14); \
    memcpy(&arg15, plasma_ptr, sizeof(arg15)); plasma_ptr += sizeof(arg15); \
    memcpy(&arg16, plasma_ptr, sizeof(arg16)); plasma_ptr += sizeof(arg16); \
}

#define plasma_unpack_args_20( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10, \
    arg11, \
    arg12, \
    arg13, \
    arg14, \
    arg15, \
    arg16, \
    arg17, \
    arg18, \
    arg19, \
    arg20) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
    memcpy(&arg11, plasma_ptr, sizeof(arg11)); plasma_ptr += sizeof(arg11); \
    memcpy(&arg12, plasma_ptr, sizeof(arg12)); plasma_ptr += sizeof(arg12); \
    memcpy(&arg13, plasma_ptr, sizeof(arg13)); plasma_ptr += sizeof(arg13); \
    memcpy(&arg14, plasma_ptr, sizeof(arg14)); plasma_ptr += sizeof(arg14); \
    memcpy(&arg15, plasma_ptr, sizeof(arg15)); plasma_ptr += sizeof(arg15); \
    memcpy(&arg16, plasma_ptr, sizeof(arg16)); plasma_ptr += sizeof(arg16); \
    memcpy(&arg17, plasma_ptr, sizeof(arg17)); plasma_ptr += sizeof(arg17); \
    memcpy(&arg18, plasma_ptr, sizeof(arg18)); plasma_ptr += sizeof(arg18); \
    memcpy(&arg19, plasma_ptr, sizeof(arg19)); plasma_ptr += sizeof(arg19); \
    memcpy(&arg20, plasma_ptr, sizeof(arg20)); plasma_ptr += sizeof(arg20); \
}

#define plasma_unpack_args_21( \
    arg1, \
    arg2, \
    arg3, \
    arg4, \
    arg5, \
    arg6, \
    arg7, \
    arg8, \
    arg9, \
    arg10, \
    arg11, \
    arg12, \
    arg13, \
    arg14, \
    arg15, \
    arg16, \
    arg17, \
    arg18, \
    arg19, \
    arg20, \
    arg21) \
{ \
    unsigned char *plasma_ptr = plasma->args_buff; \
    memcpy(&arg1, plasma_ptr, sizeof(arg1)); plasma_ptr += sizeof(arg1); \
    memcpy(&arg2, plasma_ptr, sizeof(arg2)); plasma_ptr += sizeof(arg2); \
    memcpy(&arg3, plasma_ptr, sizeof(arg3)); plasma_ptr += sizeof(arg3); \
    memcpy(&arg4, plasma_ptr, sizeof(arg4)); plasma_ptr += sizeof(arg4); \
    memcpy(&arg5, plasma_ptr, sizeof(arg5)); plasma_ptr += sizeof(arg5); \
    memcpy(&arg6, plasma_ptr, sizeof(arg6)); plasma_ptr += sizeof(arg6); \
    memcpy(&arg7, plasma_ptr, sizeof(arg7)); plasma_ptr += sizeof(arg7); \
    memcpy(&arg8, plasma_ptr, sizeof(arg8)); plasma_ptr += sizeof(arg8); \
    memcpy(&arg9, plasma_ptr, sizeof(arg9)); plasma_ptr += sizeof(arg9); \
    memcpy(&arg10, plasma_ptr, sizeof(arg10)); plasma_ptr += sizeof(arg10); \
    memcpy(&arg11, plasma_ptr, sizeof(arg11)); plasma_ptr += sizeof(arg11); \
    memcpy(&arg12, plasma_ptr, sizeof(arg12)); plasma_ptr += sizeof(arg12); \
    memcpy(&arg13, plasma_ptr, sizeof(arg13)); plasma_ptr += sizeof(arg13); \
    memcpy(&arg14, plasma_ptr, sizeof(arg14)); plasma_ptr += sizeof(arg14); \
    memcpy(&arg15, plasma_ptr, sizeof(arg15)); plasma_ptr += sizeof(arg15); \
    memcpy(&arg16, plasma_ptr, sizeof(arg16)); plasma_ptr += sizeof(arg16); \
    memcpy(&arg17, plasma_ptr, sizeof(arg17)); plasma_ptr += sizeof(arg17); \
    memcpy(&arg18, plasma_ptr, sizeof(arg18)); plasma_ptr += sizeof(arg18); \
    memcpy(&arg19, plasma_ptr, sizeof(arg19)); plasma_ptr += sizeof(arg19); \
    memcpy(&arg20, plasma_ptr, sizeof(arg20)); plasma_ptr += sizeof(arg20); \
    memcpy(&arg21, plasma_ptr, sizeof(arg21)); plasma_ptr += sizeof(arg21); \
}

#endif
