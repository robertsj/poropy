#ifndef NUCLEARDATA_H_
#define NUCLEARDATA_H_

typedef struct{
	float SM2_NO_XE;
	float M2_NO_XE;
	float M2_XE;
	float NUFISS1;
	float K2;
	float K1;
	float NUFISS2;
	float BOR1_NO_XE;
	float BOR2_XE;
	float REMOV1;
	float DIFF2;
	float DIFF1;
	float BOR1_XE;
	float XE_YIELD;
	float NU;
	float BOR2_NO_XE;
	float KAPPA;
	float ABS1;
	float ABS2;
	float XE2_MAC;
	float K_INF_NO_XE;
	float XE2_MIC;
	float SM2_XE;
	float K_INF_XE;
} NUCLEAR_DATA;

NUCLEAR_DATA get_2g_parms(float B, float ENR, int IFBA, int WABA, int GAD);

float get_SM2_NO_XE_0(float B, float R);
float get_SM2_NO_XE_1(float B, float R);
float get_M2_NO_XE_0(float B, float R);
float get_M2_NO_XE_1(float B, float R);
float get_M2_XE_0(float B, float R);
float get_M2_XE_1(float B, float R);
float get_NUFISS1_0(float B, float R);
float get_NUFISS1_1(float B, float R);
float get_K2_0(float B, float R);
float get_K2_1(float B, float R);
float get_K1_0(float B, float R);
float get_K1_1(float B, float R);
float get_NUFISS2_0(float B, float R);
float get_NUFISS2_1(float B, float R);
float get_BOR1_NO_XE_0(float B, float R);
float get_BOR1_NO_XE_1(float B, float R);
float get_BOR2_XE_0(float B, float R);
float get_BOR2_XE_1(float B, float R);
float get_REMOV1_0(float B, float R);
float get_REMOV1_1(float B, float R);
float get_DIFF2_0(float B, float R);
float get_DIFF2_1(float B, float R);
float get_DIFF1_0(float B, float R);
float get_DIFF1_1(float B, float R);
float get_BOR1_XE_0(float B, float R);
float get_BOR1_XE_1(float B, float R);
float get_XE_YIELD_0(float B, float R);
float get_XE_YIELD_1(float B, float R);
float get_NU_0(float B, float R);
float get_NU_1(float B, float R);
float get_BOR2_NO_XE_0(float B, float R);
float get_BOR2_NO_XE_1(float B, float R);
float get_KAPPA_0(float B, float R);
float get_KAPPA_1(float B, float R);
float get_ABS1_0(float B, float R);
float get_ABS1_1(float B, float R);
float get_ABS2_0(float B, float R);
float get_ABS2_1(float B, float R);
float get_XE2_MAC_0(float B, float R);
float get_XE2_MAC_1(float B, float R);
float get_K_INF_NO_XE_0(float B, float R);
float get_K_INF_NO_XE_1(float B, float R);
float get_XE2_MIC_0(float B, float R);
float get_XE2_MIC_1(float B, float R);
float get_SM2_XE_0(float B, float R);
float get_SM2_XE_1(float B, float R);
float get_K_INF_XE_0(float B, float R);
float get_K_INF_XE_1(float B, float R);

#endif /* NUCLEARDATA_H_ */