#include <math.h>
#include "nucleardata.h"

#include "nucleardata_IFBA.h"
#include "nucleardata_WABA.h"
#include "nucleardata_GAD.h"

float XE_CUTOFF = 0.1;
float IFBA_CUTOFF = 30;
float GAD_CUTOFF= 15.0;
float WABA_CUTOFF = 22.5;


NUCLEAR_DATA get_2g_parms(float B, float ENR, int IFBA, int WABA, int GAD){
    NUCLEAR_DATA r;
    if (B < XE_CUTOFF){
		r.SM2_NO_XE = get_SM2_NO_XE_0(B,ENR);
		r.M2_NO_XE = get_M2_NO_XE_0(B,ENR);
		r.M2_XE = get_M2_XE_0(B,ENR);
		r.NUFISS1 = get_NUFISS1_0(B,ENR);
		r.K2 = get_K2_0(B,ENR);
		r.K1 = get_K1_0(B,ENR);
		r.NUFISS2 = get_NUFISS2_0(B,ENR);
		r.BOR1_NO_XE = get_BOR1_NO_XE_0(B,ENR);
		r.BOR2_XE = get_BOR2_XE_0(B,ENR);
		r.REMOV1 = get_REMOV1_0(B,ENR);
		r.DIFF2 = get_DIFF2_0(B,ENR);
		r.DIFF1 = get_DIFF1_0(B,ENR);
		r.BOR1_XE = get_BOR1_XE_0(B,ENR);
		r.XE_YIELD = get_XE_YIELD_0(B,ENR);
		r.NU = get_NU_0(B,ENR);
		r.BOR2_NO_XE = get_BOR2_NO_XE_0(B,ENR);
		r.KAPPA = get_KAPPA_0(B,ENR);
		r.ABS1 = get_ABS1_0(B,ENR);
		r.ABS2 = get_ABS2_0(B,ENR);
		r.XE2_MAC = get_XE2_MAC_0(B,ENR);
		r.K_INF_NO_XE = get_K_INF_NO_XE_0(B,ENR);
		r.XE2_MIC = get_XE2_MIC_0(B,ENR);
		r.SM2_XE = get_SM2_XE_0(B,ENR);
		r.K_INF_XE = get_K_INF_XE_0(B,ENR);

    } else {
		r.SM2_NO_XE = get_SM2_NO_XE_1(B,ENR);
		r.M2_NO_XE = get_M2_NO_XE_1(B,ENR);
		r.M2_XE = get_M2_XE_1(B,ENR);
		r.NUFISS1 = get_NUFISS1_1(B,ENR);
		r.K2 = get_K2_1(B,ENR);
		r.K1 = get_K1_1(B,ENR);
		r.NUFISS2 = get_NUFISS2_1(B,ENR);
		r.BOR1_NO_XE = get_BOR1_NO_XE_1(B,ENR);
		r.BOR2_XE = get_BOR2_XE_1(B,ENR);
		r.REMOV1 = get_REMOV1_1(B,ENR);
		r.DIFF2 = get_DIFF2_1(B,ENR);
		r.DIFF1 = get_DIFF1_1(B,ENR);
		r.BOR1_XE = get_BOR1_XE_1(B,ENR);
		r.XE_YIELD = get_XE_YIELD_1(B,ENR);
		r.NU = get_NU_1(B,ENR);
		r.BOR2_NO_XE = get_BOR2_NO_XE_1(B,ENR);
		r.KAPPA = get_KAPPA_1(B,ENR);
		r.ABS1 = get_ABS1_1(B,ENR);
		r.ABS2 = get_ABS2_1(B,ENR);
		r.XE2_MAC = get_XE2_MAC_1(B,ENR);
		r.K_INF_NO_XE = get_K_INF_NO_XE_1(B,ENR);
		r.XE2_MIC = get_XE2_MIC_1(B,ENR);
		r.SM2_XE = get_SM2_XE_1(B,ENR);
		r.K_INF_XE = get_K_INF_XE_1(B,ENR);

    }
    if (IFBA) {
        if (B < XE_CUTOFF){
			r.SM2_NO_XE = r.SM2_NO_XE - get_SM2_NO_XE_IFBA_DIFF_0(B,ENR);
			r.M2_NO_XE = r.M2_NO_XE - get_M2_NO_XE_IFBA_DIFF_0(B,ENR);
			r.M2_XE = r.M2_XE - get_M2_XE_IFBA_DIFF_0(B,ENR);
			r.NUFISS1 = r.NUFISS1 - get_NUFISS1_IFBA_DIFF_0(B,ENR);
			r.K2 = r.K2 - get_K2_IFBA_DIFF_0(B,ENR);
			r.K1 = r.K1 - get_K1_IFBA_DIFF_0(B,ENR);
			r.NUFISS2 = r.NUFISS2 - get_NUFISS2_IFBA_DIFF_0(B,ENR);
			r.BOR1_NO_XE = r.BOR1_NO_XE - get_BOR1_NO_XE_IFBA_DIFF_0(B,ENR);
			r.BOR2_XE = r.BOR2_XE - get_BOR2_XE_IFBA_DIFF_0(B,ENR);
			r.REMOV1 = r.REMOV1 - get_REMOV1_IFBA_DIFF_0(B,ENR);
			r.DIFF2 = r.DIFF2 - get_DIFF2_IFBA_DIFF_0(B,ENR);
			r.DIFF1 = r.DIFF1 - get_DIFF1_IFBA_DIFF_0(B,ENR);
			r.BOR1_XE = r.BOR1_XE - get_BOR1_XE_IFBA_DIFF_0(B,ENR);
			r.XE_YIELD = r.XE_YIELD - get_XE_YIELD_IFBA_DIFF_0(B,ENR);
			r.NU = r.NU - get_NU_IFBA_DIFF_0(B,ENR);
			r.BOR2_NO_XE = r.BOR2_NO_XE - get_BOR2_NO_XE_IFBA_DIFF_0(B,ENR);
			r.KAPPA = r.KAPPA - get_KAPPA_IFBA_DIFF_0(B,ENR);
			r.ABS1 = r.ABS1 - get_ABS1_IFBA_DIFF_0(B,ENR);
			r.ABS2 = r.ABS2 - get_ABS2_IFBA_DIFF_0(B,ENR);
			r.XE2_MAC = r.XE2_MAC - get_XE2_MAC_IFBA_DIFF_0(B,ENR);
			r.K_INF_NO_XE = r.K_INF_NO_XE - get_K_INF_NO_XE_IFBA_DIFF_0(B,ENR);
			r.XE2_MIC = r.XE2_MIC - get_XE2_MIC_IFBA_DIFF_0(B,ENR);
			r.SM2_XE = r.SM2_XE - get_SM2_XE_IFBA_DIFF_0(B,ENR);
			r.K_INF_XE = r.K_INF_XE - get_K_INF_XE_IFBA_DIFF_0(B,ENR);

        } else if (B < IFBA_CUTOFF){
			r.SM2_NO_XE = r.SM2_NO_XE - get_SM2_NO_XE_IFBA_DIFF_1(B,ENR);
			r.M2_NO_XE = r.M2_NO_XE - get_M2_NO_XE_IFBA_DIFF_1(B,ENR);
			r.M2_XE = r.M2_XE - get_M2_XE_IFBA_DIFF_1(B,ENR);
			r.NUFISS1 = r.NUFISS1 - get_NUFISS1_IFBA_DIFF_1(B,ENR);
			r.K2 = r.K2 - get_K2_IFBA_DIFF_1(B,ENR);
			r.K1 = r.K1 - get_K1_IFBA_DIFF_1(B,ENR);
			r.NUFISS2 = r.NUFISS2 - get_NUFISS2_IFBA_DIFF_1(B,ENR);
			r.BOR1_NO_XE = r.BOR1_NO_XE - get_BOR1_NO_XE_IFBA_DIFF_1(B,ENR);
			r.BOR2_XE = r.BOR2_XE - get_BOR2_XE_IFBA_DIFF_1(B,ENR);
			r.REMOV1 = r.REMOV1 - get_REMOV1_IFBA_DIFF_1(B,ENR);
			r.DIFF2 = r.DIFF2 - get_DIFF2_IFBA_DIFF_1(B,ENR);
			r.DIFF1 = r.DIFF1 - get_DIFF1_IFBA_DIFF_1(B,ENR);
			r.BOR1_XE = r.BOR1_XE - get_BOR1_XE_IFBA_DIFF_1(B,ENR);
			r.XE_YIELD = r.XE_YIELD - get_XE_YIELD_IFBA_DIFF_1(B,ENR);
			r.NU = r.NU - get_NU_IFBA_DIFF_1(B,ENR);
			r.BOR2_NO_XE = r.BOR2_NO_XE - get_BOR2_NO_XE_IFBA_DIFF_1(B,ENR);
			r.KAPPA = r.KAPPA - get_KAPPA_IFBA_DIFF_1(B,ENR);
			r.ABS1 = r.ABS1 - get_ABS1_IFBA_DIFF_1(B,ENR);
			r.ABS2 = r.ABS2 - get_ABS2_IFBA_DIFF_1(B,ENR);
			r.XE2_MAC = r.XE2_MAC - get_XE2_MAC_IFBA_DIFF_1(B,ENR);
			r.K_INF_NO_XE = r.K_INF_NO_XE - get_K_INF_NO_XE_IFBA_DIFF_1(B,ENR);
			r.XE2_MIC = r.XE2_MIC - get_XE2_MIC_IFBA_DIFF_1(B,ENR);
			r.SM2_XE = r.SM2_XE - get_SM2_XE_IFBA_DIFF_1(B,ENR);
			r.K_INF_XE = r.K_INF_XE - get_K_INF_XE_IFBA_DIFF_1(B,ENR);

        } else {
			r.SM2_NO_XE = r.SM2_NO_XE - get_SM2_NO_XE_IFBA_DIFF_2(B,ENR);
			r.M2_NO_XE = r.M2_NO_XE - get_M2_NO_XE_IFBA_DIFF_2(B,ENR);
			r.M2_XE = r.M2_XE - get_M2_XE_IFBA_DIFF_2(B,ENR);
			r.NUFISS1 = r.NUFISS1 - get_NUFISS1_IFBA_DIFF_2(B,ENR);
			r.K2 = r.K2 - get_K2_IFBA_DIFF_2(B,ENR);
			r.K1 = r.K1 - get_K1_IFBA_DIFF_2(B,ENR);
			r.NUFISS2 = r.NUFISS2 - get_NUFISS2_IFBA_DIFF_2(B,ENR);
			r.BOR1_NO_XE = r.BOR1_NO_XE - get_BOR1_NO_XE_IFBA_DIFF_2(B,ENR);
			r.BOR2_XE = r.BOR2_XE - get_BOR2_XE_IFBA_DIFF_2(B,ENR);
			r.REMOV1 = r.REMOV1 - get_REMOV1_IFBA_DIFF_2(B,ENR);
			r.DIFF2 = r.DIFF2 - get_DIFF2_IFBA_DIFF_2(B,ENR);
			r.DIFF1 = r.DIFF1 - get_DIFF1_IFBA_DIFF_2(B,ENR);
			r.BOR1_XE = r.BOR1_XE - get_BOR1_XE_IFBA_DIFF_2(B,ENR);
			r.XE_YIELD = r.XE_YIELD - get_XE_YIELD_IFBA_DIFF_2(B,ENR);
			r.NU = r.NU - get_NU_IFBA_DIFF_2(B,ENR);
			r.BOR2_NO_XE = r.BOR2_NO_XE - get_BOR2_NO_XE_IFBA_DIFF_2(B,ENR);
			r.KAPPA = r.KAPPA - get_KAPPA_IFBA_DIFF_2(B,ENR);
			r.ABS1 = r.ABS1 - get_ABS1_IFBA_DIFF_2(B,ENR);
			r.ABS2 = r.ABS2 - get_ABS2_IFBA_DIFF_2(B,ENR);
			r.XE2_MAC = r.XE2_MAC - get_XE2_MAC_IFBA_DIFF_2(B,ENR);
			r.K_INF_NO_XE = r.K_INF_NO_XE - get_K_INF_NO_XE_IFBA_DIFF_2(B,ENR);
			r.XE2_MIC = r.XE2_MIC - get_XE2_MIC_IFBA_DIFF_2(B,ENR);
			r.SM2_XE = r.SM2_XE - get_SM2_XE_IFBA_DIFF_2(B,ENR);
			r.K_INF_XE = r.K_INF_XE - get_K_INF_XE_IFBA_DIFF_2(B,ENR);

        }
    }
    if (GAD){
        if (B < XE_CUTOFF){
			r.SM2_NO_XE = r.SM2_NO_XE - get_SM2_NO_XE_GAD_DIFF_0(B,ENR);
			r.M2_NO_XE = r.M2_NO_XE - get_M2_NO_XE_GAD_DIFF_0(B,ENR);
			r.M2_XE = r.M2_XE - get_M2_XE_GAD_DIFF_0(B,ENR);
			r.NUFISS1 = r.NUFISS1 - get_NUFISS1_GAD_DIFF_0(B,ENR);
			r.K2 = r.K2 - get_K2_GAD_DIFF_0(B,ENR);
			r.K1 = r.K1 - get_K1_GAD_DIFF_0(B,ENR);
			r.NUFISS2 = r.NUFISS2 - get_NUFISS2_GAD_DIFF_0(B,ENR);
			r.BOR1_NO_XE = r.BOR1_NO_XE - get_BOR1_NO_XE_GAD_DIFF_0(B,ENR);
			r.BOR2_XE = r.BOR2_XE - get_BOR2_XE_GAD_DIFF_0(B,ENR);
			r.REMOV1 = r.REMOV1 - get_REMOV1_GAD_DIFF_0(B,ENR);
			r.DIFF2 = r.DIFF2 - get_DIFF2_GAD_DIFF_0(B,ENR);
			r.DIFF1 = r.DIFF1 - get_DIFF1_GAD_DIFF_0(B,ENR);
			r.BOR1_XE = r.BOR1_XE - get_BOR1_XE_GAD_DIFF_0(B,ENR);
			r.XE_YIELD = r.XE_YIELD - get_XE_YIELD_GAD_DIFF_0(B,ENR);
			r.NU = r.NU - get_NU_GAD_DIFF_0(B,ENR);
			r.BOR2_NO_XE = r.BOR2_NO_XE - get_BOR2_NO_XE_GAD_DIFF_0(B,ENR);
			r.KAPPA = r.KAPPA - get_KAPPA_GAD_DIFF_0(B,ENR);
			r.ABS1 = r.ABS1 - get_ABS1_GAD_DIFF_0(B,ENR);
			r.ABS2 = r.ABS2 - get_ABS2_GAD_DIFF_0(B,ENR);
			r.XE2_MAC = r.XE2_MAC - get_XE2_MAC_GAD_DIFF_0(B,ENR);
			r.K_INF_NO_XE = r.K_INF_NO_XE - get_K_INF_NO_XE_GAD_DIFF_0(B,ENR);
			r.XE2_MIC = r.XE2_MIC - get_XE2_MIC_GAD_DIFF_0(B,ENR);
			r.SM2_XE = r.SM2_XE - get_SM2_XE_GAD_DIFF_0(B,ENR);
			r.K_INF_XE = r.K_INF_XE - get_K_INF_XE_GAD_DIFF_0(B,ENR);

        } else if (B < GAD_CUTOFF){
			r.SM2_NO_XE = r.SM2_NO_XE - get_SM2_NO_XE_GAD_DIFF_1(B,ENR);
			r.M2_NO_XE = r.M2_NO_XE - get_M2_NO_XE_GAD_DIFF_1(B,ENR);
			r.M2_XE = r.M2_XE - get_M2_XE_GAD_DIFF_1(B,ENR);
			r.NUFISS1 = r.NUFISS1 - get_NUFISS1_GAD_DIFF_1(B,ENR);
			r.K2 = r.K2 - get_K2_GAD_DIFF_1(B,ENR);
			r.K1 = r.K1 - get_K1_GAD_DIFF_1(B,ENR);
			r.NUFISS2 = r.NUFISS2 - get_NUFISS2_GAD_DIFF_1(B,ENR);
			r.BOR1_NO_XE = r.BOR1_NO_XE - get_BOR1_NO_XE_GAD_DIFF_1(B,ENR);
			r.BOR2_XE = r.BOR2_XE - get_BOR2_XE_GAD_DIFF_1(B,ENR);
			r.REMOV1 = r.REMOV1 - get_REMOV1_GAD_DIFF_1(B,ENR);
			r.DIFF2 = r.DIFF2 - get_DIFF2_GAD_DIFF_1(B,ENR);
			r.DIFF1 = r.DIFF1 - get_DIFF1_GAD_DIFF_1(B,ENR);
			r.BOR1_XE = r.BOR1_XE - get_BOR1_XE_GAD_DIFF_1(B,ENR);
			r.XE_YIELD = r.XE_YIELD - get_XE_YIELD_GAD_DIFF_1(B,ENR);
			r.NU = r.NU - get_NU_GAD_DIFF_1(B,ENR);
			r.BOR2_NO_XE = r.BOR2_NO_XE - get_BOR2_NO_XE_GAD_DIFF_1(B,ENR);
			r.KAPPA = r.KAPPA - get_KAPPA_GAD_DIFF_1(B,ENR);
			r.ABS1 = r.ABS1 - get_ABS1_GAD_DIFF_1(B,ENR);
			r.ABS2 = r.ABS2 - get_ABS2_GAD_DIFF_1(B,ENR);
			r.XE2_MAC = r.XE2_MAC - get_XE2_MAC_GAD_DIFF_1(B,ENR);
			r.K_INF_NO_XE = r.K_INF_NO_XE - get_K_INF_NO_XE_GAD_DIFF_1(B,ENR);
			r.XE2_MIC = r.XE2_MIC - get_XE2_MIC_GAD_DIFF_1(B,ENR);
			r.SM2_XE = r.SM2_XE - get_SM2_XE_GAD_DIFF_1(B,ENR);
			r.K_INF_XE = r.K_INF_XE - get_K_INF_XE_GAD_DIFF_1(B,ENR);

        } else {
			r.SM2_NO_XE = r.SM2_NO_XE - get_SM2_NO_XE_GAD_DIFF_2(B,ENR);
			r.M2_NO_XE = r.M2_NO_XE - get_M2_NO_XE_GAD_DIFF_2(B,ENR);
			r.M2_XE = r.M2_XE - get_M2_XE_GAD_DIFF_2(B,ENR);
			r.NUFISS1 = r.NUFISS1 - get_NUFISS1_GAD_DIFF_2(B,ENR);
			r.K2 = r.K2 - get_K2_GAD_DIFF_2(B,ENR);
			r.K1 = r.K1 - get_K1_GAD_DIFF_2(B,ENR);
			r.NUFISS2 = r.NUFISS2 - get_NUFISS2_GAD_DIFF_2(B,ENR);
			r.BOR1_NO_XE = r.BOR1_NO_XE - get_BOR1_NO_XE_GAD_DIFF_2(B,ENR);
			r.BOR2_XE = r.BOR2_XE - get_BOR2_XE_GAD_DIFF_2(B,ENR);
			r.REMOV1 = r.REMOV1 - get_REMOV1_GAD_DIFF_2(B,ENR);
			r.DIFF2 = r.DIFF2 - get_DIFF2_GAD_DIFF_2(B,ENR);
			r.DIFF1 = r.DIFF1 - get_DIFF1_GAD_DIFF_2(B,ENR);
			r.BOR1_XE = r.BOR1_XE - get_BOR1_XE_GAD_DIFF_2(B,ENR);
			r.XE_YIELD = r.XE_YIELD - get_XE_YIELD_GAD_DIFF_2(B,ENR);
			r.NU = r.NU - get_NU_GAD_DIFF_2(B,ENR);
			r.BOR2_NO_XE = r.BOR2_NO_XE - get_BOR2_NO_XE_GAD_DIFF_2(B,ENR);
			r.KAPPA = r.KAPPA - get_KAPPA_GAD_DIFF_2(B,ENR);
			r.ABS1 = r.ABS1 - get_ABS1_GAD_DIFF_2(B,ENR);
			r.ABS2 = r.ABS2 - get_ABS2_GAD_DIFF_2(B,ENR);
			r.XE2_MAC = r.XE2_MAC - get_XE2_MAC_GAD_DIFF_2(B,ENR);
			r.K_INF_NO_XE = r.K_INF_NO_XE - get_K_INF_NO_XE_GAD_DIFF_2(B,ENR);
			r.XE2_MIC = r.XE2_MIC - get_XE2_MIC_GAD_DIFF_2(B,ENR);
			r.SM2_XE = r.SM2_XE - get_SM2_XE_GAD_DIFF_2(B,ENR);
			r.K_INF_XE = r.K_INF_XE - get_K_INF_XE_GAD_DIFF_2(B,ENR);

        }
    }
    if (WABA){
        if (B < XE_CUTOFF){
			r.SM2_NO_XE = r.SM2_NO_XE - get_SM2_NO_XE_WABA_DIFF_0(B,ENR);
			r.M2_NO_XE = r.M2_NO_XE - get_M2_NO_XE_WABA_DIFF_0(B,ENR);
			r.M2_XE = r.M2_XE - get_M2_XE_WABA_DIFF_0(B,ENR);
			r.NUFISS1 = r.NUFISS1 - get_NUFISS1_WABA_DIFF_0(B,ENR);
			r.K2 = r.K2 - get_K2_WABA_DIFF_0(B,ENR);
			r.K1 = r.K1 - get_K1_WABA_DIFF_0(B,ENR);
			r.NUFISS2 = r.NUFISS2 - get_NUFISS2_WABA_DIFF_0(B,ENR);
			r.BOR1_NO_XE = r.BOR1_NO_XE - get_BOR1_NO_XE_WABA_DIFF_0(B,ENR);
			r.BOR2_XE = r.BOR2_XE - get_BOR2_XE_WABA_DIFF_0(B,ENR);
			r.REMOV1 = r.REMOV1 - get_REMOV1_WABA_DIFF_0(B,ENR);
			r.DIFF2 = r.DIFF2 - get_DIFF2_WABA_DIFF_0(B,ENR);
			r.DIFF1 = r.DIFF1 - get_DIFF1_WABA_DIFF_0(B,ENR);
			r.BOR1_XE = r.BOR1_XE - get_BOR1_XE_WABA_DIFF_0(B,ENR);
			r.XE_YIELD = r.XE_YIELD - get_XE_YIELD_WABA_DIFF_0(B,ENR);
			r.NU = r.NU - get_NU_WABA_DIFF_0(B,ENR);
			r.BOR2_NO_XE = r.BOR2_NO_XE - get_BOR2_NO_XE_WABA_DIFF_0(B,ENR);
			r.KAPPA = r.KAPPA - get_KAPPA_WABA_DIFF_0(B,ENR);
			r.ABS1 = r.ABS1 - get_ABS1_WABA_DIFF_0(B,ENR);
			r.ABS2 = r.ABS2 - get_ABS2_WABA_DIFF_0(B,ENR);
			r.XE2_MAC = r.XE2_MAC - get_XE2_MAC_WABA_DIFF_0(B,ENR);
			r.K_INF_NO_XE = r.K_INF_NO_XE - get_K_INF_NO_XE_WABA_DIFF_0(B,ENR);
			r.XE2_MIC = r.XE2_MIC - get_XE2_MIC_WABA_DIFF_0(B,ENR);
			r.SM2_XE = r.SM2_XE - get_SM2_XE_WABA_DIFF_0(B,ENR);
			r.K_INF_XE = r.K_INF_XE - get_K_INF_XE_WABA_DIFF_0(B,ENR);

        } else if (B < WABA_CUTOFF){
			r.SM2_NO_XE = r.SM2_NO_XE - get_SM2_NO_XE_WABA_DIFF_1(B,ENR);
			r.M2_NO_XE = r.M2_NO_XE - get_M2_NO_XE_WABA_DIFF_1(B,ENR);
			r.M2_XE = r.M2_XE - get_M2_XE_WABA_DIFF_1(B,ENR);
			r.NUFISS1 = r.NUFISS1 - get_NUFISS1_WABA_DIFF_1(B,ENR);
			r.K2 = r.K2 - get_K2_WABA_DIFF_1(B,ENR);
			r.K1 = r.K1 - get_K1_WABA_DIFF_1(B,ENR);
			r.NUFISS2 = r.NUFISS2 - get_NUFISS2_WABA_DIFF_1(B,ENR);
			r.BOR1_NO_XE = r.BOR1_NO_XE - get_BOR1_NO_XE_WABA_DIFF_1(B,ENR);
			r.BOR2_XE = r.BOR2_XE - get_BOR2_XE_WABA_DIFF_1(B,ENR);
			r.REMOV1 = r.REMOV1 - get_REMOV1_WABA_DIFF_1(B,ENR);
			r.DIFF2 = r.DIFF2 - get_DIFF2_WABA_DIFF_1(B,ENR);
			r.DIFF1 = r.DIFF1 - get_DIFF1_WABA_DIFF_1(B,ENR);
			r.BOR1_XE = r.BOR1_XE - get_BOR1_XE_WABA_DIFF_1(B,ENR);
			r.XE_YIELD = r.XE_YIELD - get_XE_YIELD_WABA_DIFF_1(B,ENR);
			r.NU = r.NU - get_NU_WABA_DIFF_1(B,ENR);
			r.BOR2_NO_XE = r.BOR2_NO_XE - get_BOR2_NO_XE_WABA_DIFF_1(B,ENR);
			r.KAPPA = r.KAPPA - get_KAPPA_WABA_DIFF_1(B,ENR);
			r.ABS1 = r.ABS1 - get_ABS1_WABA_DIFF_1(B,ENR);
			r.ABS2 = r.ABS2 - get_ABS2_WABA_DIFF_1(B,ENR);
			r.XE2_MAC = r.XE2_MAC - get_XE2_MAC_WABA_DIFF_1(B,ENR);
			r.K_INF_NO_XE = r.K_INF_NO_XE - get_K_INF_NO_XE_WABA_DIFF_1(B,ENR);
			r.XE2_MIC = r.XE2_MIC - get_XE2_MIC_WABA_DIFF_1(B,ENR);
			r.SM2_XE = r.SM2_XE - get_SM2_XE_WABA_DIFF_1(B,ENR);
			r.K_INF_XE = r.K_INF_XE - get_K_INF_XE_WABA_DIFF_1(B,ENR);

        } else {
			r.SM2_NO_XE = r.SM2_NO_XE - get_SM2_NO_XE_WABA_DIFF_2(B,ENR);
			r.M2_NO_XE = r.M2_NO_XE - get_M2_NO_XE_WABA_DIFF_2(B,ENR);
			r.M2_XE = r.M2_XE - get_M2_XE_WABA_DIFF_2(B,ENR);
			r.NUFISS1 = r.NUFISS1 - get_NUFISS1_WABA_DIFF_2(B,ENR);
			r.K2 = r.K2 - get_K2_WABA_DIFF_2(B,ENR);
			r.K1 = r.K1 - get_K1_WABA_DIFF_2(B,ENR);
			r.NUFISS2 = r.NUFISS2 - get_NUFISS2_WABA_DIFF_2(B,ENR);
			r.BOR1_NO_XE = r.BOR1_NO_XE - get_BOR1_NO_XE_WABA_DIFF_2(B,ENR);
			r.BOR2_XE = r.BOR2_XE - get_BOR2_XE_WABA_DIFF_2(B,ENR);
			r.REMOV1 = r.REMOV1 - get_REMOV1_WABA_DIFF_2(B,ENR);
			r.DIFF2 = r.DIFF2 - get_DIFF2_WABA_DIFF_2(B,ENR);
			r.DIFF1 = r.DIFF1 - get_DIFF1_WABA_DIFF_2(B,ENR);
			r.BOR1_XE = r.BOR1_XE - get_BOR1_XE_WABA_DIFF_2(B,ENR);
			r.XE_YIELD = r.XE_YIELD - get_XE_YIELD_WABA_DIFF_2(B,ENR);
			r.NU = r.NU - get_NU_WABA_DIFF_2(B,ENR);
			r.BOR2_NO_XE = r.BOR2_NO_XE - get_BOR2_NO_XE_WABA_DIFF_2(B,ENR);
			r.KAPPA = r.KAPPA - get_KAPPA_WABA_DIFF_2(B,ENR);
			r.ABS1 = r.ABS1 - get_ABS1_WABA_DIFF_2(B,ENR);
			r.ABS2 = r.ABS2 - get_ABS2_WABA_DIFF_2(B,ENR);
			r.XE2_MAC = r.XE2_MAC - get_XE2_MAC_WABA_DIFF_2(B,ENR);
			r.K_INF_NO_XE = r.K_INF_NO_XE - get_K_INF_NO_XE_WABA_DIFF_2(B,ENR);
			r.XE2_MIC = r.XE2_MIC - get_XE2_MIC_WABA_DIFF_2(B,ENR);
			r.SM2_XE = r.SM2_XE - get_SM2_XE_WABA_DIFF_2(B,ENR);
			r.K_INF_XE = r.K_INF_XE - get_K_INF_XE_WABA_DIFF_2(B,ENR);

        }
    }

    return r;
}

//SM2 NO XE
float c0[] = {-4.49810650009e-07, 8.83425672229e-06, -6.65638881746e-05, 0.000198676652178, 0.000664380292657, 1.45684307382e-24, -2.77286777595e-23, 2.43790494817e-22, -1.71083460889e-21, 1.6558052386e-20};
float get_SM2_NO_XE_0(float B, float R){
	float val = (c0[0]*pow(R,4) + c0[1]*pow(R,3) + c0[2]*pow(R,2) + c0[3]*R + c0[4])*B + (c0[5]*pow(R,4) + c0[6]*pow(R,3) + c0[7]*pow(R,2) + c0[8]*R + c0[9]);
	return val;
}
//SM2 NO XE
float c1[] = {-3.76928530072e-13, 7.33276281017e-12, -6.44526827371e-11, 1.53285532914e-10, -3.52405422008e-10, 4.54204118976e-11, -8.41552502978e-10, 7.49991821636e-09, -1.48283977742e-08, 4.83407178199e-08, -2.42243791126e-09, 4.12353831043e-08, -3.3415497721e-07, 4.24396417432e-07, -2.36285317177e-06, 6.6750911908e-08, -1.16474316502e-06, 8.77563044815e-06, -1.22957898652e-05, 6.31440662551e-05, -2.05191608598e-07, 3.9820255855e-06, -3.63869544429e-05, 0.000259706261794, -5.15983919492e-05};
float get_SM2_NO_XE_1(float B, float R){
	float val = (c1[0]*pow(R,4) + c1[1]*pow(R,3) + c1[2]*pow(R,2) + c1[3]*R + c1[4])*pow(B,4) + (c1[5]*pow(R,4) + c1[6]*pow(R,3) + c1[7]*pow(R,2) + c1[8]*R + c1[9])*pow(B,3) + (c1[10]*pow(R,4) + c1[11]*pow(R,3) + c1[12]*pow(R,2) + c1[13]*R + c1[14])*pow(B,2) + (c1[15]*pow(R,4) + c1[16]*pow(R,3) + c1[17]*pow(R,2) + c1[18]*R + c1[19])*B + (c1[20]*pow(R,4) + c1[21]*pow(R,3) + c1[22]*pow(R,2) + c1[23]*R + c1[24]);
	return val;
}
//M2 NO XE
float c2[] = {0.0182111370132, -0.298662582369, 1.88303014006, -5.28757448368, 1.5218622181, -0.00109260345321, 0.0256886763589, -0.237975053179, 0.761107045172, 60.9675917517};
float get_M2_NO_XE_0(float B, float R){
	float val = (c2[0]*pow(R,4) + c2[1]*pow(R,3) + c2[2]*pow(R,2) + c2[3]*R + c2[4])*B + (c2[5]*pow(R,4) + c2[6]*pow(R,3) + c2[7]*pow(R,2) + c2[8]*R + c2[9]);
	return val;
}
//M2 NO XE
float c3[] = {3.4714273055e-09, -5.86245002062e-08, 3.70713302908e-07, -1.1166671541e-06, 1.68031677916e-06, -5.9031478239e-07, 1.0161491481e-05, -6.58820814505e-05, 0.000203256024383, -0.00030858875762, 3.12341151242e-05, -0.000544754825567, 0.00361819389785, -0.0116590197721, 0.019290646342, -0.000481785906119, 0.00830393464773, -0.055271423312, 0.193774003624, -0.462067901024, 0.000348161622287, 0.00038438016966, -0.070258735219, 0.283881796852, 61.0099766412};
float get_M2_NO_XE_1(float B, float R){
	float val = (c3[0]*pow(R,4) + c3[1]*pow(R,3) + c3[2]*pow(R,2) + c3[3]*R + c3[4])*pow(B,4) + (c3[5]*pow(R,4) + c3[6]*pow(R,3) + c3[7]*pow(R,2) + c3[8]*R + c3[9])*pow(B,3) + (c3[10]*pow(R,4) + c3[11]*pow(R,3) + c3[12]*pow(R,2) + c3[13]*R + c3[14])*pow(B,2) + (c3[15]*pow(R,4) + c3[16]*pow(R,3) + c3[17]*pow(R,2) + c3[18]*R + c3[19])*B + (c3[20]*pow(R,4) + c3[21]*pow(R,3) + c3[22]*pow(R,2) + c3[23]*R + c3[24]);
	return val;
}
//M2 XE
float c4[] = {0.0236741400303, -0.398454054885, 2.52555824103, -6.79842947588, 0.843610688663, -0.00109260345321, 0.0256886763589, -0.237975053179, 0.761107045172, 60.9675917517};
float get_M2_XE_0(float B, float R){
	float val = (c4[0]*pow(R,4) + c4[1]*pow(R,3) + c4[2]*pow(R,2) + c4[3]*R + c4[4])*B + (c4[5]*pow(R,4) + c4[6]*pow(R,3) + c4[7]*pow(R,2) + c4[8]*R + c4[9]);
	return val;
}
//M2 XE
float c5[] = {4.06911920561e-09, -6.77012035605e-08, 4.20824061775e-07, -1.23007384584e-06, 1.74896404504e-06, -6.32207376949e-07, 1.07616073465e-05, -6.89094871793e-05, 0.000208641792508, -0.000307425542979, 2.99701234476e-05, -0.000520937538175, 0.00344972313677, -0.0110908530876, 0.0183912507293, -0.000400264276184, 0.006925042046, -0.0464818279012, 0.167783906672, -0.428117796396, 0.000146322533086, 0.00388319560935, -0.0946371918065, 0.384242879538, 60.6850530004};
float get_M2_XE_1(float B, float R){
	float val = (c5[0]*pow(R,4) + c5[1]*pow(R,3) + c5[2]*pow(R,2) + c5[3]*R + c5[4])*pow(B,4) + (c5[5]*pow(R,4) + c5[6]*pow(R,3) + c5[7]*pow(R,2) + c5[8]*R + c5[9])*pow(B,3) + (c5[10]*pow(R,4) + c5[11]*pow(R,3) + c5[12]*pow(R,2) + c5[13]*R + c5[14])*pow(B,2) + (c5[15]*pow(R,4) + c5[16]*pow(R,3) + c5[17]*pow(R,2) + c5[18]*R + c5[19])*B + (c5[20]*pow(R,4) + c5[21]*pow(R,3) + c5[22]*pow(R,2) + c5[23]*R + c5[24]);
	return val;
}
//NUFISS1
float c6[] = {-5.28115686208e-06, 8.83592328681e-05, -0.000553924359729, 0.00158819720384, -0.0017240860653, -8.48624277579e-06, 0.00013899957825, -0.000866015286859, 0.00359587826743, 0.000380133406295};
float get_NUFISS1_0(float B, float R){
	float val = (c6[0]*pow(R,4) + c6[1]*pow(R,3) + c6[2]*pow(R,2) + c6[3]*R + c6[4])*B + (c6[5]*pow(R,4) + c6[6]*pow(R,3) + c6[7]*pow(R,2) + c6[8]*R + c6[9]);
	return val;
}
//NUFISS1
float c7[] = {-5.21232210374e-09, 9.85274081292e-08, -6.91803304329e-07, 2.07293614222e-06, -1.78582989102e-06, 3.18817430786e-07, -6.3054394911e-06, 4.76452151641e-05, -0.000166982814727, 0.000152608046711, -1.67462066933e-06, 3.69520754632e-05, -0.000318829350046, 0.00238450462037, 0.00128583870001};
float get_NUFISS1_1(float B, float R){
	float val = (c7[0]*pow(R,4) + c7[1]*pow(R,3) + c7[2]*pow(R,2) + c7[3]*R + c7[4])*pow(B,2) + (c7[5]*pow(R,4) + c7[6]*pow(R,3) + c7[7]*pow(R,2) + c7[8]*R + c7[9])*B + (c7[10]*pow(R,4) + c7[11]*pow(R,3) + c7[12]*pow(R,2) + c7[13]*R + c7[14]);
	return val;
}
//K2
float c8[] = {0.000455266301703, -0.00914177792205, 0.0729880562325, -0.258552003393, -0.0455654830311, -0.000209415606237, 0.00568480840817, -0.0588501912184, 0.262769200352, 0.554620384853};
float get_K2_0(float B, float R){
	float val = (c8[0]*pow(R,4) + c8[1]*pow(R,3) + c8[2]*pow(R,2) + c8[3]*R + c8[4])*B + (c8[5]*pow(R,4) + c8[6]*pow(R,3) + c8[7]*pow(R,2) + c8[8]*R + c8[9]);
	return val;
}
//K2
float c9[] = {-3.05372821352e-07, 5.17029838595e-06, -3.00821341311e-05, 5.93191869535e-05, 2.18599290901e-05, 2.75694282809e-05, -0.00050982016344, 0.00345891233262, -0.0092639802761, 0.000551662636378, -0.000384579307195, 0.00808434326162, -0.0687776067618, 0.269815579333, 0.532848926372};
float get_K2_1(float B, float R){
	float val = (c9[0]*pow(R,4) + c9[1]*pow(R,3) + c9[2]*pow(R,2) + c9[3]*R + c9[4])*pow(B,2) + (c9[5]*pow(R,4) + c9[6]*pow(R,3) + c9[7]*pow(R,2) + c9[8]*R + c9[9])*B + (c9[10]*pow(R,4) + c9[11]*pow(R,3) + c9[12]*pow(R,2) + c9[13]*R + c9[14]);
	return val;
}
//K1
float c10[] = {-0.000127482887564, 0.00213927952437, -0.0131132064986, 0.0344486084537, -0.0482439183713, -0.000338725728989, 0.00557404064546, -0.035361928837, 0.149434403827, 0.0025313454318};
float get_K1_0(float B, float R){
	float val = (c10[0]*pow(R,4) + c10[1]*pow(R,3) + c10[2]*pow(R,2) + c10[3]*R + c10[4])*B + (c10[5]*pow(R,4) + c10[6]*pow(R,3) + c10[7]*pow(R,2) + c10[8]*R + c10[9]);
	return val;
}
//K1
float c11[] = {-2.08779823145e-07, 3.9517167771e-06, -2.78811455836e-05, 8.4453850636e-05, -7.19533393143e-05, 1.29051773846e-05, -0.000255040371851, 0.00193607152447, -0.00685814232696, 0.00612476558632, -6.76625163706e-05, 0.00149048634847, -0.0132053456125, 0.0989647441614, 0.0405435674915};
float get_K1_1(float B, float R){
	float val = (c11[0]*pow(R,4) + c11[1]*pow(R,3) + c11[2]*pow(R,2) + c11[3]*R + c11[4])*pow(B,2) + (c11[5]*pow(R,4) + c11[6]*pow(R,3) + c11[7]*pow(R,2) + c11[8]*R + c11[9])*B + (c11[10]*pow(R,4) + c11[11]*pow(R,3) + c11[12]*pow(R,2) + c11[13]*R + c11[14]);
	return val;
}
//NUFISS2
float c12[] = {-0.000127476210006, 0.00212460370812, -0.0129260894173, 0.0286829564356, -0.0309604838766, -1.09288421435e-05, 0.000264702066003, -0.00349569773068, 0.0503164270374, 0.00177188706315};
float get_NUFISS2_0(float B, float R){
	float val = (c12[0]*pow(R,4) + c12[1]*pow(R,3) + c12[2]*pow(R,2) + c12[3]*R + c12[4])*B + (c12[5]*pow(R,4) + c12[6]*pow(R,3) + c12[7]*pow(R,2) + c12[8]*R + c12[9]);
	return val;
}
//NUFISS2
float c13[] = {4.74290910538e-11, -7.37176859244e-10, 3.60875319587e-09, -3.00158462414e-09, -1.7941278799e-08, -7.28948186031e-09, 1.21383821337e-07, -6.8184207714e-07, 1.08280268901e-06, 2.11763581737e-06, 2.68956219875e-07, -4.67105282989e-06, 2.84289297524e-05, -5.47459274777e-05, -0.000104113748093, -6.75254779857e-07, 3.0312847545e-06, 8.04233698123e-05, -0.000944169919526, 0.00418086143604, 1.40257754928e-06, 4.93340105271e-05, -0.00201583588249, 0.0448102878772, 0.00910199008969};
float get_NUFISS2_1(float B, float R){
	float val = (c13[0]*pow(R,4) + c13[1]*pow(R,3) + c13[2]*pow(R,2) + c13[3]*R + c13[4])*pow(B,4) + (c13[5]*pow(R,4) + c13[6]*pow(R,3) + c13[7]*pow(R,2) + c13[8]*R + c13[9])*pow(B,3) + (c13[10]*pow(R,4) + c13[11]*pow(R,3) + c13[12]*pow(R,2) + c13[13]*R + c13[14])*pow(B,2) + (c13[15]*pow(R,4) + c13[16]*pow(R,3) + c13[17]*pow(R,2) + c13[18]*R + c13[19])*B + (c13[20]*pow(R,4) + c13[21]*pow(R,3) + c13[22]*pow(R,2) + c13[23]*R + c13[24]);
	return val;
}
//BOR1 NO XE
float c14[] = {1.09266174051e-07, -1.64142102012e-06, 8.97876947693e-06, -2.18219740765e-05, 3.34213078762e-05, 5.46337578789e-08, -9.91175867129e-07, 7.15524195172e-06, -3.02238611059e-05, 0.000245939910787};
float get_BOR1_NO_XE_0(float B, float R){
	float val = (c14[0]*pow(R,4) + c14[1]*pow(R,3) + c14[2]*pow(R,2) + c14[3]*R + c14[4])*B + (c14[5]*pow(R,4) + c14[6]*pow(R,3) + c14[7]*pow(R,2) + c14[8]*R + c14[9]);
	return val;
}
//BOR1 NO XE
float c15[] = {1.01047259144e-13, -1.71626815723e-12, 1.08732251095e-11, -3.07835915262e-11, 3.29839767286e-11, -1.52993055214e-11, 2.59300064136e-10, -1.63865216179e-09, 4.63595377612e-09, -5.01520021242e-09, 7.8882208031e-10, -1.33824401636e-08, 8.46066058908e-08, -2.39572235281e-07, 2.62475416827e-07, -1.52554039501e-08, 2.63394243928e-07, -1.70550425324e-06, 5.01077877601e-06, -5.85464694193e-06, 6.46105400061e-08, -1.19287258907e-06, 8.60237355306e-06, -3.47227862624e-05, 0.000252815513691};
float get_BOR1_NO_XE_1(float B, float R){
	float val = (c15[0]*pow(R,4) + c15[1]*pow(R,3) + c15[2]*pow(R,2) + c15[3]*R + c15[4])*pow(B,4) + (c15[5]*pow(R,4) + c15[6]*pow(R,3) + c15[7]*pow(R,2) + c15[8]*R + c15[9])*pow(B,3) + (c15[10]*pow(R,4) + c15[11]*pow(R,3) + c15[12]*pow(R,2) + c15[13]*R + c15[14])*pow(B,2) + (c15[15]*pow(R,4) + c15[16]*pow(R,3) + c15[17]*pow(R,2) + c15[18]*R + c15[19])*B + (c15[20]*pow(R,4) + c15[21]*pow(R,3) + c15[22]*pow(R,2) + c15[23]*R + c15[24]);
	return val;
}
//BOR2 XE
float c16[] = {-0.0728426337013, 1.00522591818, -4.74639281369, 8.31847873065, -14.8633626414, 0.0163848029432, -0.319217795815, 2.89892169653, -16.9159833494, 452.119085912};
float get_BOR2_XE_0(float B, float R){
	float val = (c16[0]*pow(R,4) + c16[1]*pow(R,3) + c16[2]*pow(R,2) + c16[3]*R + c16[4])*B + (c16[5]*pow(R,4) + c16[6]*pow(R,3) + c16[7]*pow(R,2) + c16[8]*R + c16[9]);
	return val;
}
//BOR2 XE
float c17[] = {-1.04646810948e-08, 1.71399082417e-07, -9.59794345469e-07, 2.02753841129e-06, -1.35477162464e-06, 1.38005891418e-06, -2.38269234861e-05, 0.000143541160545, -0.000341774129782, 0.000289801213355, -2.95729900999e-05, 0.000550500691661, -0.00351361165227, 0.00896532764234, -0.0114982341691, -0.00109217730495, 0.0208538983001, -0.161059832232, 0.579510327628, -0.424271423897, 0.00866704689804, -0.205857982877, 2.32083823996, -15.7909467448, 450.235530622};
float get_BOR2_XE_1(float B, float R){
	float val = (c17[0]*pow(R,4) + c17[1]*pow(R,3) + c17[2]*pow(R,2) + c17[3]*R + c17[4])*pow(B,4) + (c17[5]*pow(R,4) + c17[6]*pow(R,3) + c17[7]*pow(R,2) + c17[8]*R + c17[9])*pow(B,3) + (c17[10]*pow(R,4) + c17[11]*pow(R,3) + c17[12]*pow(R,2) + c17[13]*R + c17[14])*pow(B,2) + (c17[15]*pow(R,4) + c17[16]*pow(R,3) + c17[17]*pow(R,2) + c17[18]*R + c17[19])*B + (c17[20]*pow(R,4) + c17[21]*pow(R,3) + c17[22]*pow(R,2) + c17[23]*R + c17[24]);
	return val;
}
//REMOV1
float c18[] = {3.27795907001e-05, -0.000521074137193, 0.00306277985797, -0.00796013035459, 0.00893845797355, 5.09923183417e-06, -9.25628502113e-05, 0.000669858318466, -0.00285625696867, 0.0205351912513};
float get_REMOV1_0(float B, float R){
	float val = (c18[0]*pow(R,4) + c18[1]*pow(R,3) + c18[2]*pow(R,2) + c18[3]*R + c18[4])*B + (c18[5]*pow(R,4) + c18[6]*pow(R,3) + c18[7]*pow(R,2) + c18[8]*R + c18[9]);
	return val;
}
//REMOV1
float c19[] = {9.40657765025e-12, -1.60753904865e-10, 1.0291564301e-09, -2.9665259572e-09, 3.25830081169e-09, -1.47884666396e-09, 2.5180061705e-08, -1.60639542501e-07, 4.62688990367e-07, -5.14679525361e-07, 7.83272645417e-08, -1.33319491182e-06, 8.50332151224e-06, -2.45477435932e-05, 2.79229562382e-05, -1.51675973445e-06, 2.62077013141e-05, -0.000170729639555, 0.000511042239808, -0.000634902402688, 6.25500150533e-06, -0.000114720660828, 0.00082296045283, -0.00332194901242, 0.0212409051783};
float get_REMOV1_1(float B, float R){
	float val = (c19[0]*pow(R,4) + c19[1]*pow(R,3) + c19[2]*pow(R,2) + c19[3]*R + c19[4])*pow(B,4) + (c19[5]*pow(R,4) + c19[6]*pow(R,3) + c19[7]*pow(R,2) + c19[8]*R + c19[9])*pow(B,3) + (c19[10]*pow(R,4) + c19[11]*pow(R,3) + c19[12]*pow(R,2) + c19[13]*R + c19[14])*pow(B,2) + (c19[15]*pow(R,4) + c19[16]*pow(R,3) + c19[17]*pow(R,2) + c19[18]*R + c19[19])*B + (c19[20]*pow(R,4) + c19[21]*pow(R,3) + c19[22]*pow(R,2) + c19[23]*R + c19[24]);
	return val;
}
//DIFF2
float c20[] = {7.28449868776e-05, -0.00103440091678, 0.00522738340717, -0.0110958601666, 0.0120639163837, -1.27510832387e-05, 0.000253308340666, -0.00217784501246, 0.0077035964039, 0.372064027347};
float get_DIFF2_0(float B, float R){
	float val = (c20[0]*pow(R,4) + c20[1]*pow(R,3) + c20[2]*pow(R,2) + c20[3]*R + c20[4])*B + (c20[5]*pow(R,4) + c20[6]*pow(R,3) + c20[7]*pow(R,2) + c20[8]*R + c20[9]);
	return val;
}
//DIFF2
float c21[] = {2.68328863608e-06, -3.92624734243e-05, 3.63164340157e-05, -0.000213951387211, 0.369092881568, 1.08224570273e-07, -1.99889103116e-06, 1.5781290783e-05, -4.23638711631e-05, -1.96774949102e-05, -1.15355797539e-05, 0.000236576519403, -0.00193482200693, 0.00736894179019, 0.00364937297609, -0.000100082722932, 0.00217902690573, -0.0186978105104, 0.0813944500507, -0.225769818424};
float get_DIFF2_1(float B, float R){
	float val = (c21[0]*pow(R,4) + c21[1]*pow(R,3) + c21[2]*pow(R,2) + c21[3]*R + c21[4])*exp(B*(c21[5]*pow(R,4) + c21[6]*pow(R,3) + c21[7]*pow(R,2) + c21[8]*R + c21[9])) + (c21[10]*pow(R,4) + c21[11]*pow(R,3) + c21[12]*pow(R,2) + c21[13]*R + c21[14])*exp(B*(c21[15]*pow(R,4) + c21[16]*pow(R,3) + c21[17]*pow(R,2) + c21[18]*R + c21[19]));
	return val;
}
//DIFF1
float c22[] = {0.000546388907325, -0.00905787999315, 0.0558696611103, -0.151284406713, 0.131408458841, -1.82150913032e-05, 0.000449266216344, -0.00430032124641, 0.0233430291955, 1.39922043742};
float get_DIFF1_0(float B, float R){
	float val = (c22[0]*pow(R,4) + c22[1]*pow(R,3) + c22[2]*pow(R,2) + c22[3]*R + c22[4])*B + (c22[5]*pow(R,4) + c22[6]*pow(R,3) + c22[7]*pow(R,2) + c22[8]*R + c22[9]);
	return val;
}
//DIFF1
float c23[] = {4.68951327282e-11, -6.52055424002e-10, 3.13417883613e-09, -5.59161472598e-09, 2.11726059161e-09, -6.38242037836e-09, 8.76151256819e-08, -4.13053937421e-07, 7.00524477507e-07, -1.84689000076e-07, 2.22352078499e-07, -2.76610719888e-06, 1.07043802594e-05, -8.36002410479e-06, -1.53660285964e-05, -7.76652147193e-07, -8.09202037424e-06, 0.000205957725021, -0.00111770093408, 0.00191507181414, 8.5677315868e-06, 3.32955627415e-05, -0.00187296504056, 0.017027295364, 1.40332443961};
float get_DIFF1_1(float B, float R){
	float val = (c23[0]*pow(R,4) + c23[1]*pow(R,3) + c23[2]*pow(R,2) + c23[3]*R + c23[4])*pow(B,4) + (c23[5]*pow(R,4) + c23[6]*pow(R,3) + c23[7]*pow(R,2) + c23[8]*R + c23[9])*pow(B,3) + (c23[10]*pow(R,4) + c23[11]*pow(R,3) + c23[12]*pow(R,2) + c23[13]*R + c23[14])*pow(B,2) + (c23[15]*pow(R,4) + c23[16]*pow(R,3) + c23[17]*pow(R,2) + c23[18]*R + c23[19])*B + (c23[20]*pow(R,4) + c23[21]*pow(R,3) + c23[22]*pow(R,2) + c23[23]*R + c23[24]);
	return val;
}
//BOR1 XE
float c24[] = {0.000546230070659, -0.00784127029818, 0.0350503206536, -0.0775870666695, 0.715582159686, 0.00298659774567, -0.0534673825547, 0.378537397204, -1.5385929285, 11.7822457239};
float get_BOR1_XE_0(float B, float R){
	float val = (c24[0]*pow(R,4) + c24[1]*pow(R,3) + c24[2]*pow(R,2) + c24[3]*R + c24[4])*B + (c24[5]*pow(R,4) + c24[6]*pow(R,3) + c24[7]*pow(R,2) + c24[8]*R + c24[9]);
	return val;
}
//BOR1 XE
float c25[] = {4.46908167637e-09, -7.62839064342e-08, 4.8567722434e-07, -1.38241340925e-06, 1.48963938288e-06, -6.8909936786e-07, 1.17179549177e-05, -7.43011180885e-05, 0.000211016449147, -0.000229313126216, 3.61583955364e-05, -0.000614667894928, 0.00389418747417, -0.0110542091171, 0.0121509063933, -0.000711648925444, 0.0123035464294, -0.0797816920897, 0.234808743377, -0.274981832835, 0.00319697854568, -0.058803569683, 0.421855148867, -1.68507048173, 12.0425084327};
float get_BOR1_XE_1(float B, float R){
	float val = (c25[0]*pow(R,4) + c25[1]*pow(R,3) + c25[2]*pow(R,2) + c25[3]*R + c25[4])*pow(B,4) + (c25[5]*pow(R,4) + c25[6]*pow(R,3) + c25[7]*pow(R,2) + c25[8]*R + c25[9])*pow(B,3) + (c25[10]*pow(R,4) + c25[11]*pow(R,3) + c25[12]*pow(R,2) + c25[13]*R + c25[14])*pow(B,2) + (c25[15]*pow(R,4) + c25[16]*pow(R,3) + c25[17]*pow(R,2) + c25[18]*R + c25[19])*B + (c25[20]*pow(R,4) + c25[21]*pow(R,3) + c25[22]*pow(R,2) + c25[23]*R + c25[24]);
	return val;
}
//XE-YIELD
float c26[] = {1.09271767908e-05, -0.000174349205052, 0.00104682307615, -0.00284807124364, 0.00315408438073, -1.27481108898e-06, 1.95471964718e-05, -0.000108870058616, 0.000247833539475, 0.0654778642278};
float get_XE_YIELD_0(float B, float R){
	float val = (c26[0]*pow(R,4) + c26[1]*pow(R,3) + c26[2]*pow(R,2) + c26[3]*R + c26[4])*B + (c26[5]*pow(R,4) + c26[6]*pow(R,3) + c26[7]*pow(R,2) + c26[8]*R + c26[9]);
	return val;
}
//XE-YIELD
float c27[] = {-1.60372721757e-13, 1.97789310996e-12, -5.21393985802e-12, -1.96979445334e-11, 8.97778497187e-11, 2.12088472059e-11, -2.26906191793e-10, 1.0671488777e-10, 6.12762671683e-09, -1.99573437781e-08, -4.93528699638e-10, -9.37873661919e-10, 9.93860121403e-08, -7.38862739024e-07, 1.78725743834e-06, -2.11522173873e-08, 6.84749202997e-07, -7.61545014165e-06, 3.7845810806e-05, -7.80498721994e-05, 4.68611427437e-07, -1.25391002666e-05, 0.000128581932188, -0.00063809636072, 0.00149499374857, 2.36081037131e-06, -4.28054942233e-05, 0.000295122856686, -0.000933951419569, 0.0668069911302};
float get_XE_YIELD_1(float B, float R){
	float val = (c27[0]*pow(R,4) + c27[1]*pow(R,3) + c27[2]*pow(R,2) + c27[3]*R + c27[4])*pow(B,5) + (c27[5]*pow(R,4) + c27[6]*pow(R,3) + c27[7]*pow(R,2) + c27[8]*R + c27[9])*pow(B,4) + (c27[10]*pow(R,4) + c27[11]*pow(R,3) + c27[12]*pow(R,2) + c27[13]*R + c27[14])*pow(B,3) + (c27[15]*pow(R,4) + c27[16]*pow(R,3) + c27[17]*pow(R,2) + c27[18]*R + c27[19])*pow(B,2) + (c27[20]*pow(R,4) + c27[21]*pow(R,3) + c27[22]*pow(R,2) + c27[23]*R + c27[24])*B + (c27[25]*pow(R,4) + c27[26]*pow(R,3) + c27[27]*pow(R,2) + c27[28]*R + c27[29]);
	return val;
}
//NU
float c28[] = {-0.00236738388085, 0.0377567325122, -0.223065765681, 0.576748772887, -0.53929314031, 0.000109269357726, -0.00182115451413, 0.0114208104824, -0.0323948845462, 2.49581911256};
float get_NU_0(float B, float R){
	float val = (c28[0]*pow(R,4) + c28[1]*pow(R,3) + c28[2]*pow(R,2) + c28[3]*R + c28[4])*B + (c28[5]*pow(R,4) + c28[6]*pow(R,3) + c28[7]*pow(R,2) + c28[8]*R + c28[9]);
	return val;
}
//NU
float c29[] = {-7.54904007328e-12, 1.07250831441e-10, -4.53762205474e-10, 1.44360374036e-10, 2.22658421525e-09, 1.0504407836e-09, -1.42594303292e-08, 5.07368404734e-08, 6.30955822409e-08, -5.22183861387e-07, -2.67220098106e-08, 1.9503346182e-07, 1.57841073031e-06, -1.91009812227e-05, 5.37954340856e-05, -1.06270885023e-06, 2.86213612931e-05, -0.000292031169824, 0.00138788012767, -0.00278976855125, 2.80901670408e-05, -0.000653898345135, 0.00609427166181, -0.0282586741348, 0.0629520687515, 4.66192133207e-05, -0.000987981478793, 0.00781359372811, -0.0279493731345, 2.49910967805};
float get_NU_1(float B, float R){
	float val = (c29[0]*pow(R,4) + c29[1]*pow(R,3) + c29[2]*pow(R,2) + c29[3]*R + c29[4])*pow(B,5) + (c29[5]*pow(R,4) + c29[6]*pow(R,3) + c29[7]*pow(R,2) + c29[8]*R + c29[9])*pow(B,4) + (c29[10]*pow(R,4) + c29[11]*pow(R,3) + c29[12]*pow(R,2) + c29[13]*R + c29[14])*pow(B,3) + (c29[15]*pow(R,4) + c29[16]*pow(R,3) + c29[17]*pow(R,2) + c29[18]*R + c29[19])*pow(B,2) + (c29[20]*pow(R,4) + c29[21]*pow(R,3) + c29[22]*pow(R,2) + c29[23]*R + c29[24])*B + (c29[25]*pow(R,4) + c29[26]*pow(R,3) + c29[27]*pow(R,2) + c29[28]*R + c29[29]);
	return val;
}
//BOR2 NO XE
float c30[] = {1.81927801363e-07, -4.44053750462e-06, 4.14669571644e-05, -0.000168074232905, -2.9544524325e-06, 1.82116306721e-07, -4.32217646726e-06, 4.80295072742e-05, -0.000325436096168, 0.00949329566252};
float get_BOR2_NO_XE_0(float B, float R){
	float val = (c30[0]*pow(R,4) + c30[1]*pow(R,3) + c30[2]*pow(R,2) + c30[3]*R + c30[4])*B + (c30[5]*pow(R,4) + c30[6]*pow(R,3) + c30[7]*pow(R,2) + c30[8]*R + c30[9]);
	return val;
}
//BOR2 NO XE
float c31[] = {-7.79589614344e-14, 1.18461494075e-12, -4.92572682786e-12, 4.02598646748e-13, 1.47672803904e-11, 9.62778664147e-12, -1.69458709272e-10, 9.19340567095e-10, -1.35500606647e-09, 1.00952021602e-10, 1.60032641742e-10, -1.97795369026e-09, 1.29742068538e-08, -5.54045965574e-08, 1.11575832658e-08, -3.0112547518e-08, 5.7208133133e-07, -4.29690584009e-06, 1.48754972386e-05, -1.18230780581e-05, 1.51959429755e-07, -3.92332973297e-06, 4.68739212014e-05, -0.000328431001364, 0.00947748022095};
float get_BOR2_NO_XE_1(float B, float R){
	float val = (c31[0]*pow(R,4) + c31[1]*pow(R,3) + c31[2]*pow(R,2) + c31[3]*R + c31[4])*pow(B,4) + (c31[5]*pow(R,4) + c31[6]*pow(R,3) + c31[7]*pow(R,2) + c31[8]*R + c31[9])*pow(B,3) + (c31[10]*pow(R,4) + c31[11]*pow(R,3) + c31[12]*pow(R,2) + c31[13]*R + c31[14])*pow(B,2) + (c31[15]*pow(R,4) + c31[16]*pow(R,3) + c31[17]*pow(R,2) + c31[18]*R + c31[19])*B + (c31[20]*pow(R,4) + c31[21]*pow(R,3) + c31[22]*pow(R,2) + c31[23]*R + c31[24]);
	return val;
}
//KAPPA
float c32[] = {9.10570820506e-15, -1.45934156779e-13, 8.69413066724e-13, -2.28635341269e-12, 2.26763190825e-12, 5.46006961162e-16, -8.71184365586e-15, 5.21622009682e-14, -1.40924173034e-13, 3.26032569422e-11};
float get_KAPPA_0(float B, float R){
	float val = (c32[0]*pow(R,4) + c32[1]*pow(R,3) + c32[2]*pow(R,2) + c32[3]*R + c32[4])*B + (c32[5]*pow(R,4) + c32[6]*pow(R,3) + c32[7]*pow(R,2) + c32[8]*R + c32[9]);
	return val;
}
//KAPPA
float c33[] = {2.2044574769e-18, -3.32922395458e-17, 1.64749626722e-16, -2.05568932321e-16, -4.09581963564e-16, -1.59229706364e-16, 2.55803957118e-15, -1.41959330947e-14, 2.48353888686e-14, 3.2787082089e-14, 1.67532669987e-15, -2.99950062323e-14, 2.05583963461e-13, -6.5246980553e-13, 3.33109959097e-11};
float get_KAPPA_1(float B, float R){
	float val = (c33[0]*pow(R,4) + c33[1]*pow(R,3) + c33[2]*pow(R,2) + c33[3]*R + c33[4])*pow(B,2) + (c33[5]*pow(R,4) + c33[6]*pow(R,3) + c33[7]*pow(R,2) + c33[8]*R + c33[9])*B + (c33[10]*pow(R,4) + c33[11]*pow(R,3) + c33[12]*pow(R,2) + c33[13]*R + c33[14]);
	return val;
}
//ABS1
float c34[] = {-4.00610977729e-06, 6.67686543196e-05, -0.000417001847617, 0.00117205451093, -0.000857349671919, -3.84249527629e-06, 6.21573739453e-05, -0.000375031247445, 0.0015696379314, 0.00657839472069};
float get_ABS1_0(float B, float R){
	float val = (c34[0]*pow(R,4) + c34[1]*pow(R,3) + c34[2]*pow(R,2) + c34[3]*R + c34[4])*B + (c34[5]*pow(R,4) + c34[6]*pow(R,3) + c34[7]*pow(R,2) + c34[8]*R + c34[9]);
	return val;
}
//ABS1
float c35[] = {-2.37865623753e-09, 4.47107437071e-08, -3.28213270679e-07, 1.17451354236e-06, -2.08560777104e-06, 1.29409964692e-07, -2.519894603e-06, 1.95317552069e-05, -8.04442079717e-05, 0.00020166141774, -6.17953200499e-07, 1.13410573807e-05, -7.86825344468e-05, 0.000794713913275, 0.00740213956063};
float get_ABS1_1(float B, float R){
	float val = (c35[0]*pow(R,4) + c35[1]*pow(R,3) + c35[2]*pow(R,2) + c35[3]*R + c35[4])*pow(B,2) + (c35[5]*pow(R,4) + c35[6]*pow(R,3) + c35[7]*pow(R,2) + c35[8]*R + c35[9])*B + (c35[10]*pow(R,4) + c35[11]*pow(R,3) + c35[12]*pow(R,2) + c35[13]*R + c35[14]);
	return val;
}
//ABS2
float c36[] = {0.000136588291976, -0.00207056242734, 0.0110003654624, -0.0208652123974, 0.0344754211405, -6.00810787307e-06, 0.000135878135845, -0.00165605536297, 0.02325865804, 0.0267476672046};
float get_ABS2_0(float B, float R){
	float val = (c36[0]*pow(R,4) + c36[1]*pow(R,3) + c36[2]*pow(R,2) + c36[3]*R + c36[4])*B + (c36[5]*pow(R,4) + c36[6]*pow(R,3) + c36[7]*pow(R,2) + c36[8]*R + c36[9]);
	return val;
}
//ABS2
float c37[] = {6.91033722625e-06, -0.00033879507338, 0.00354770153048, 0.00700258542531, 0.0611454119202, -2.51986393542e-06, 1.21158297627e-05, 0.0004236990971, -0.00449586549271, 0.00840770574311, -1.42754403645e-06, 0.000283622152648, -0.004093081466, 0.0141469354321, -0.0314075643878, 0.00102586576004, -0.0168896187334, 0.0955161242255, -0.179204951639, -0.100128850589};
float get_ABS2_1(float B, float R){
	float val = (c37[0]*pow(R,4) + c37[1]*pow(R,3) + c37[2]*pow(R,2) + c37[3]*R + c37[4])*exp(B*(c37[5]*pow(R,4) + c37[6]*pow(R,3) + c37[7]*pow(R,2) + c37[8]*R + c37[9])) + (c37[10]*pow(R,4) + c37[11]*pow(R,3) + c37[12]*pow(R,2) + c37[13]*R + c37[14])*exp(B*(c37[15]*pow(R,4) + c37[16]*pow(R,3) + c37[17]*pow(R,2) + c37[18]*R + c37[19]));
	return val;
}
//XE2 MAC
float c38[] = {-1.18369021591e-05, 0.000247738341948, -0.00252884726066, 0.0167631772781, -0.000881947727264, 2.18508503431e-23, -5.57946952825e-22, 6.72996846455e-21, -5.89187783253e-20, 5.50370358268e-19};
float get_XE2_MAC_0(float B, float R){
	float val = (c38[0]*pow(R,4) + c38[1]*pow(R,3) + c38[2]*pow(R,2) + c38[3]*R + c38[4])*B + (c38[5]*pow(R,4) + c38[6]*pow(R,3) + c38[7]*pow(R,2) + c38[8]*R + c38[9]);
	return val;
}
//XE2 MAC
float c39[] = {8.43446921529e-06, -0.000100075181334, 0.000445582653817, 5.73758969654e-06, 0.00176034956779, 2.10108419257e-05, -0.000458678106265, 0.00375356446818, -0.0144756397621, 0.0159148271812, -9.8549881497e-06, 0.000128667337238, -0.000724112651816, 0.00177027818778, -0.00201312588378, -0.000202464760905, 0.00696273288218, -0.0819502485605, 0.425362570077, -0.872730567624};
float get_XE2_MAC_1(float B, float R){
	float val = (c39[0]*pow(R,4) + c39[1]*pow(R,3) + c39[2]*pow(R,2) + c39[3]*R + c39[4])*exp(B*(c39[5]*pow(R,4) + c39[6]*pow(R,3) + c39[7]*pow(R,2) + c39[8]*R + c39[9])) + (c39[10]*pow(R,4) + c39[11]*pow(R,3) + c39[12]*pow(R,2) + c39[13]*R + c39[14])*exp(B*(c39[15]*pow(R,4) + c39[16]*pow(R,3) + c39[17]*pow(R,2) + c39[18]*R + c39[19]));
	return val;
}
//K-INF NO XE
float c40[] = {-5.46297205248e-05, 0.000910497384286, -0.0061134824582, 0.0233678067677, -0.0739537865285, -0.000546326876652, 0.0112329744251, -0.0940697599613, 0.41184132024, 0.557512550155};
float get_K_INF_NO_XE_0(float B, float R){
	float val = (c40[0]*pow(R,4) + c40[1]*pow(R,3) + c40[2]*pow(R,2) + c40[3]*R + c40[4])*B + (c40[5]*pow(R,4) + c40[6]*pow(R,3) + c40[7]*pow(R,2) + c40[8]*R + c40[9]);
	return val;
}
//K-INF NO XE
float c41[] = {-5.37013678191e-07, 9.53450733827e-06, -6.0572956817e-05, 0.000150016546173, -5.19361200396e-05, 4.26293452125e-05, -0.000806439888012, 0.00568957649383, -0.0169681328793, 0.00708148923869, -0.000486103944877, 0.0102936858528, -0.0880807748341, 0.391804234713, 0.577096510952};
float get_K_INF_NO_XE_1(float B, float R){
	float val = (c41[0]*pow(R,4) + c41[1]*pow(R,3) + c41[2]*pow(R,2) + c41[3]*R + c41[4])*pow(B,2) + (c41[5]*pow(R,4) + c41[6]*pow(R,3) + c41[7]*pow(R,2) + c41[8]*R + c41[9])*B + (c41[10]*pow(R,4) + c41[11]*pow(R,3) + c41[12]*pow(R,2) + c41[13]*R + c41[14]);
	return val;
}
//XE2 MIC
float c42[] = {364.252048087, -7139.22145523, 55320.7263072, -205334.50221, 83371.5370377, 127.494522066, -2668.79748403, 26791.4682913, -205893.025804, 1831789.13786};
float get_XE2_MIC_0(float B, float R){
	float val = (c42[0]*pow(R,4) + c42[1]*pow(R,3) + c42[2]*pow(R,2) + c42[3]*R + c42[4])*B + (c42[5]*pow(R,4) + c42[6]*pow(R,3) + c42[7]*pow(R,2) + c42[8]*R + c42[9]);
	return val;
}
//XE2 MIC
float c43[] = {-0.000196325402899, 0.00301305131267, -0.0156754243209, 0.0280964548748, 0.00206907553195, 0.0305177784436, -0.48784069279, 2.73263993824, -5.78549234823, 1.52212778024, -1.17886921344, 18.9676466577, -108.330741338, 241.109300293, -63.0578318965, 3.98361449442, -10.3532423788, -513.017166507, 4098.28846202, -7815.31090635, 22.7226028574, -1040.70243926, 17947.5254105, -188050.39974, 1799235.83293};
float get_XE2_MIC_1(float B, float R){
	float val = (c43[0]*pow(R,4) + c43[1]*pow(R,3) + c43[2]*pow(R,2) + c43[3]*R + c43[4])*pow(B,4) + (c43[5]*pow(R,4) + c43[6]*pow(R,3) + c43[7]*pow(R,2) + c43[8]*R + c43[9])*pow(B,3) + (c43[10]*pow(R,4) + c43[11]*pow(R,3) + c43[12]*pow(R,2) + c43[13]*R + c43[14])*pow(B,2) + (c43[15]*pow(R,4) + c43[16]*pow(R,3) + c43[17]*pow(R,2) + c43[18]*R + c43[19])*B + (c43[20]*pow(R,4) + c43[21]*pow(R,3) + c43[22]*pow(R,2) + c43[23]*R + c43[24]);
	return val;
}
//SM2 XE
float c44[] = {0.00207717395854, -18.0014014128, 357.132365551, -2434.82098702, -824.347979378, 3.46010497147, -70.4889073791, 683.051390581, -5345.64167485, 54516.4294715};
float get_SM2_XE_0(float B, float R){
	float val = (c44[0]*pow(R,4) + c44[1]*pow(R,3) + c44[2]*pow(R,2) + c44[3]*R + c44[4])*B + (c44[5]*pow(R,4) + c44[6]*pow(R,3) + c44[7]*pow(R,2) + c44[8]*R + c44[9]);
	return val;
}
//SM2 XE
float c45[] = {-7.01953975485e-07, 5.58242650782e-06, 4.35500940604e-05, -0.000509858739432, 0.00139556659239, 0.000119098391242, -0.00154217589956, 0.0019432092905, 0.0401419922501, -0.16933410465, 0.00199566397653, -0.0470192970111, 0.543760165241, -3.07397489103, 8.94739856128, -0.402828933723, 8.10841853607, -65.5811748992, 254.973603168, -390.205323799, 1.55907335701, -40.6092140421, 524.653222921, -5086.94023581, 53838.0127726};
float get_SM2_XE_1(float B, float R){
	float val = (c45[0]*pow(R,4) + c45[1]*pow(R,3) + c45[2]*pow(R,2) + c45[3]*R + c45[4])*pow(B,4) + (c45[5]*pow(R,4) + c45[6]*pow(R,3) + c45[7]*pow(R,2) + c45[8]*R + c45[9])*pow(B,3) + (c45[10]*pow(R,4) + c45[11]*pow(R,3) + c45[12]*pow(R,2) + c45[13]*R + c45[14])*pow(B,2) + (c45[15]*pow(R,4) + c45[16]*pow(R,3) + c45[17]*pow(R,2) + c45[18]*R + c45[19])*B + (c45[20]*pow(R,4) + c45[21]*pow(R,3) + c45[22]*pow(R,2) + c45[23]*R + c45[24]);
	return val;
}
//K-INF XE
float c46[] = {0.00030958745968, -0.00676235422969, 0.0587202485268, -0.22174903536, -0.0954421365789, -0.000546326876652, 0.0112329744251, -0.0940697599613, 0.41184132024, 0.557512550155};
float get_K_INF_XE_0(float B, float R){
	float val = (c46[0]*pow(R,4) + c46[1]*pow(R,3) + c46[2]*pow(R,2) + c46[3]*R + c46[4])*B + (c46[5]*pow(R,4) + c46[6]*pow(R,3) + c46[7]*pow(R,2) + c46[8]*R + c46[9]);
	return val;
}
//K-INF XE
float c47[] = {-5.01113246943e-07, 8.90885190171e-06, -5.66675803644e-05, 0.000140304234737, -4.66438938206e-05, 3.97201275447e-05, -0.000752500190387, 0.00531972391273, -0.0159203970918, 0.00647569248671, -0.000450759234143, 0.00954882312368, -0.0818157321749, 0.368313178677, 0.573870878092};
float get_K_INF_XE_1(float B, float R){
	float val = (c47[0]*pow(R,4) + c47[1]*pow(R,3) + c47[2]*pow(R,2) + c47[3]*R + c47[4])*pow(B,2) + (c47[5]*pow(R,4) + c47[6]*pow(R,3) + c47[7]*pow(R,2) + c47[8]*R + c47[9])*B + (c47[10]*pow(R,4) + c47[11]*pow(R,3) + c47[12]*pow(R,2) + c47[13]*R + c47[14]);
	return val;
}
