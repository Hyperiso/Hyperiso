#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "marty/core/looptools_interface.h"
#include "c_C7.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_C7(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t G_F = param->G_F;
    const creal_t M_W = param->M_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_cb = param->V_cb;
    const creal_t V_tb = param->V_tb;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_ub = param->V_ub;
    const ccomplex_t IT_0000 = pow(G_F, -1);
    const ccomplex_t IT_0001 = pow(m_b, -1);
    const ccomplex_t IT_0002 = pow(V_tb, -1);
    const ccomplex_t IT_0003 = cpow(conj(V_ts), -1);
    const ccomplex_t IT_0004 = pow(e_em, -1);
    const ccomplex_t IT_0005 = 9.86960440108936*IT_0000*IT_0001*IT_0002
      *IT_0003*IT_0004;
    const ccomplex_t IT_0006 = 0.101321183642338*m_t;
    const ccomplex_t IT_0007 = pow(M_W, -1);
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = pow(m_b, 2);
    const ccomplex_t IT_0020 = pow(m_s, 2);
    const ccomplex_t IT_0021 = pow(m_t, 2);
    const ccomplex_t IT_0022 = pow(m_Hp, 2);
    const ccomplex_t IT_0023 = mtylt_C0iC(0, IT_0019, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0020, IT_0021, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (-2)*IT_0024;
    const ccomplex_t IT_0028 = IT_0023*IT_0027;
    const ccomplex_t IT_0029 = mtylt_C0iC(3, IT_0019, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0020, IT_0021, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = -IT_0028 + -IT_0030;
    const ccomplex_t IT_0032 = IT_0026 + IT_0031;
    const ccomplex_t IT_0033 = IT_0018*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_u*V_ub*e_em
      *IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = 0.101321183642338*m_u;
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = pow(m_u, 2);
    const ccomplex_t IT_0042 = mtylt_C0iC(0, IT_0019, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0020, IT_0041, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0043 = IT_0025*IT_0042;
    const ccomplex_t IT_0044 = IT_0027*IT_0042;
    const ccomplex_t IT_0045 = mtylt_C0iC(3, IT_0019, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0020, IT_0041, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0046 = IT_0027*IT_0045;
    const ccomplex_t IT_0047 = -IT_0044 + -IT_0046;
    const ccomplex_t IT_0048 = IT_0043 + IT_0047;
    const ccomplex_t IT_0049 = IT_0040*IT_0048;
    const ccomplex_t IT_0050 = pow(m_c, 2);
    const ccomplex_t IT_0051 = mtylt_C0iC(0, IT_0019, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0020, IT_0050, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0052 = IT_0025*IT_0051;
    const ccomplex_t IT_0053 = IT_0027*IT_0051;
    const ccomplex_t IT_0054 = mtylt_C0iC(3, IT_0019, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0020, IT_0050, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0055 = IT_0027*IT_0054;
    const ccomplex_t IT_0056 = -IT_0053 + -IT_0055;
    const ccomplex_t IT_0057 = IT_0052 + IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = 0.101321183642338*m_c;
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = IT_0057*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = IT_0006*IT_0070;
    const ccomplex_t IT_0072 = IT_0032*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_b*V_ub*e_em
      *IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = IT_0039*IT_0077;
    const ccomplex_t IT_0079 = IT_0048*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0007*IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = IT_0063*IT_0084;
    const ccomplex_t IT_0086 = IT_0057*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*m_u*V_ub*e_em
      *IT_0007*IT_0012;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0007*IT_0012;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = IT_0039*IT_0091;
    const ccomplex_t IT_0093 = pow(M_W, 2);
    const ccomplex_t IT_0094 = mtylt_C0iC(0, IT_0020, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0019, IT_0041, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0095 = IT_0025*IT_0094;
    const ccomplex_t IT_0096 = IT_0027*IT_0094;
    const ccomplex_t IT_0097 = mtylt_C0iC(6, IT_0020, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0019, IT_0041, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0098 = IT_0027*IT_0097;
    const ccomplex_t IT_0099 = -IT_0096 + -IT_0098;
    const ccomplex_t IT_0100 = IT_0095 + IT_0099;
    const ccomplex_t IT_0101 = IT_0092*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0007*IT_0012;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0007*IT_0012;
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = IT_0063*IT_0106;
    const ccomplex_t IT_0108 = mtylt_C0iC(0, IT_0020, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0019, IT_0050, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0109 = IT_0025*IT_0108;
    const ccomplex_t IT_0110 = IT_0027*IT_0108;
    const ccomplex_t IT_0111 = mtylt_C0iC(6, IT_0020, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0019, IT_0050, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0112 = IT_0027*IT_0111;
    const ccomplex_t IT_0113 = -IT_0110 + -IT_0112;
    const ccomplex_t IT_0114 = IT_0109 + IT_0113;
    const ccomplex_t IT_0115 = IT_0107*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0007*IT_0012;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0007*IT_0012;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = IT_0117*IT_0119;
    const ccomplex_t IT_0121 = IT_0006*IT_0120;
    const ccomplex_t IT_0122 = mtylt_C0iC(0, IT_0020, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0019, IT_0021, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0123 = IT_0025*IT_0122;
    const ccomplex_t IT_0124 = IT_0027*IT_0122;
    const ccomplex_t IT_0125 = mtylt_C0iC(6, IT_0020, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0019, IT_0021, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0126 = IT_0027*IT_0125;
    const ccomplex_t IT_0127 = -IT_0124 + -IT_0126;
    const ccomplex_t IT_0128 = IT_0123 + IT_0127;
    const ccomplex_t IT_0129 = IT_0121*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0007*IT_0012;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1.4142135623731)*m_b*V_ub*e_em
      *IT_0007*IT_0012;
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = IT_0131*IT_0133;
    const ccomplex_t IT_0135 = IT_0039*IT_0134;
    const ccomplex_t IT_0136 = IT_0100*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0007*IT_0012;
    const ccomplex_t IT_0138 = 0.5*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0007*IT_0012;
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = IT_0063*IT_0141;
    const ccomplex_t IT_0143 = IT_0114*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0007*IT_0012;
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0007*IT_0012;
    const ccomplex_t IT_0147 = (-0.5)*IT_0146;
    const ccomplex_t IT_0148 = IT_0145*IT_0147;
    const ccomplex_t IT_0149 = IT_0006*IT_0148;
    const ccomplex_t IT_0150 = IT_0128*IT_0149;
    const ccomplex_t IT_0151 = IT_0033 + IT_0049 + IT_0065 + IT_0072 + IT_0079
       + IT_0086 + IT_0101 + IT_0115 + IT_0129 + IT_0136 + IT_0143 + IT_0150;
    const ccomplex_t IT_0152 = mtylt_C0iC(0, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0041, IT_0041, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0153 = mtylt_C0iC(3, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0041, IT_0041, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0154 = mtylt_C0iC(6, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0041, IT_0041, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0155 = IT_0152 + IT_0153 + IT_0154;
    const ccomplex_t IT_0156 = 0.666666666666667*IT_0024;
    const ccomplex_t IT_0157 = IT_0074*IT_0076*IT_0156;
    const ccomplex_t IT_0158 = IT_0039*IT_0157;
    const ccomplex_t IT_0159 = IT_0155*IT_0158;
    const ccomplex_t IT_0160 = mtylt_C0iC(3, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0041, IT_0093, IT_0041, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0161 = mtylt_C0iC(6, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0041, IT_0093, IT_0041, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0162 = mtylt_C0iC(0, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0041, IT_0093, IT_0041, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0163 = IT_0160 + IT_0161 + IT_0162;
    const ccomplex_t IT_0164 = IT_0131*IT_0133*IT_0156;
    const ccomplex_t IT_0165 = IT_0039*IT_0164;
    const ccomplex_t IT_0166 = IT_0163*IT_0165;
    const ccomplex_t IT_0167 = IT_0059*IT_0061*IT_0156;
    const ccomplex_t IT_0168 = IT_0063*IT_0167;
    const ccomplex_t IT_0169 = mtylt_C0iC(0, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0050, IT_0050, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0170 = mtylt_C0iC(3, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0050, IT_0050, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0171 = mtylt_C0iC(6, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0050, IT_0050, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0172 = IT_0169 + IT_0170 + IT_0171;
    const ccomplex_t IT_0173 = IT_0168*IT_0172;
    const ccomplex_t IT_0174 = mtylt_C0iC(3, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0050, IT_0093, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0175 = mtylt_C0iC(6, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0050, IT_0093, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0176 = mtylt_C0iC(0, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0050, IT_0093, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0177 = IT_0174 + IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = IT_0103*IT_0105*IT_0156;
    const ccomplex_t IT_0179 = IT_0063*IT_0178;
    const ccomplex_t IT_0180 = IT_0177*IT_0179;
    const ccomplex_t IT_0181 = IT_0081*IT_0083*IT_0156;
    const ccomplex_t IT_0182 = IT_0063*IT_0181;
    const ccomplex_t IT_0183 = IT_0172*IT_0182;
    const ccomplex_t IT_0184 = IT_0138*IT_0140*IT_0156;
    const ccomplex_t IT_0185 = IT_0063*IT_0184;
    const ccomplex_t IT_0186 = IT_0177*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*M_W*e_em;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1.4142135623731)*V_ub*e_em
      *IT_0012;
    const ccomplex_t IT_0189 = 0.5*IT_0188;
    const ccomplex_t IT_0190 = IT_0090*IT_0187*IT_0189;
    const ccomplex_t IT_0191 = 0.101321183642338*IT_0190;
    const ccomplex_t IT_0192 = IT_0097*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0012;
    const ccomplex_t IT_0194 = 0.5*IT_0193;
    const ccomplex_t IT_0195 = IT_0105*IT_0187*IT_0194;
    const ccomplex_t IT_0196 = 0.101321183642338*IT_0195;
    const ccomplex_t IT_0197 = IT_0111*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0012;
    const ccomplex_t IT_0199 = 0.5*IT_0198;
    const ccomplex_t IT_0200 = IT_0119*IT_0187*IT_0199;
    const ccomplex_t IT_0201 = 0.101321183642338*IT_0200;
    const ccomplex_t IT_0202 = IT_0125*IT_0201;
    const ccomplex_t IT_0203 = IT_0035*IT_0074;
    const ccomplex_t IT_0204 = 0.101321183642338*IT_0203;
    const ccomplex_t IT_0205 = mtylt_C0iC(6, IT_0019, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0020, IT_0041, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0206 = IT_0025*IT_0205;
    const ccomplex_t IT_0207 = m_s*IT_0206;
    const ccomplex_t IT_0208 = IT_0025*IT_0045;
    const ccomplex_t IT_0209 = m_b*IT_0208;
    const ccomplex_t IT_0210 = IT_0207 + IT_0209;
    const ccomplex_t IT_0211 = IT_0027*IT_0205;
    const ccomplex_t IT_0212 = m_s*IT_0211;
    const ccomplex_t IT_0213 = mtylt_C0iC(15, IT_0019, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0020, IT_0041, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0214 = IT_0027*IT_0213;
    const ccomplex_t IT_0215 = m_s*IT_0214;
    const ccomplex_t IT_0216 = m_b*IT_0046;
    const ccomplex_t IT_0217 = mtylt_C0iC(12, IT_0019, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0020, IT_0041, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0218 = IT_0027*IT_0217;
    const ccomplex_t IT_0219 = m_b*IT_0218;
    const ccomplex_t IT_0220 = -IT_0212 + -IT_0215 + -IT_0216 + -IT_0219;
    const ccomplex_t IT_0221 = IT_0210 + IT_0220;
    const ccomplex_t IT_0222 = IT_0204*IT_0221;
    const ccomplex_t IT_0223 = IT_0059*IT_0081;
    const ccomplex_t IT_0224 = 0.101321183642338*IT_0223;
    const ccomplex_t IT_0225 = mtylt_C0iC(6, IT_0019, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0020, IT_0050, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0226 = IT_0025*IT_0225;
    const ccomplex_t IT_0227 = m_s*IT_0226;
    const ccomplex_t IT_0228 = IT_0025*IT_0054;
    const ccomplex_t IT_0229 = m_b*IT_0228;
    const ccomplex_t IT_0230 = IT_0227 + IT_0229;
    const ccomplex_t IT_0231 = IT_0027*IT_0225;
    const ccomplex_t IT_0232 = m_s*IT_0231;
    const ccomplex_t IT_0233 = mtylt_C0iC(15, IT_0019, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0020, IT_0050, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0234 = IT_0027*IT_0233;
    const ccomplex_t IT_0235 = m_s*IT_0234;
    const ccomplex_t IT_0236 = m_b*IT_0055;
    const ccomplex_t IT_0237 = mtylt_C0iC(12, IT_0019, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0020, IT_0050, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0238 = IT_0027*IT_0237;
    const ccomplex_t IT_0239 = m_b*IT_0238;
    const ccomplex_t IT_0240 = -IT_0232 + -IT_0235 + -IT_0236 + -IT_0239;
    const ccomplex_t IT_0241 = IT_0230 + IT_0240;
    const ccomplex_t IT_0242 = IT_0224*IT_0241;
    const ccomplex_t IT_0243 = IT_0014*IT_0067;
    const ccomplex_t IT_0244 = 0.101321183642338*IT_0243;
    const ccomplex_t IT_0245 = mtylt_C0iC(6, IT_0019, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0020, IT_0021, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0246 = IT_0025*IT_0245;
    const ccomplex_t IT_0247 = m_s*IT_0246;
    const ccomplex_t IT_0248 = IT_0025*IT_0029;
    const ccomplex_t IT_0249 = m_b*IT_0248;
    const ccomplex_t IT_0250 = IT_0247 + IT_0249;
    const ccomplex_t IT_0251 = IT_0027*IT_0245;
    const ccomplex_t IT_0252 = m_s*IT_0251;
    const ccomplex_t IT_0253 = mtylt_C0iC(12, IT_0019, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0020, IT_0021, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0254 = IT_0027*IT_0253;
    const ccomplex_t IT_0255 = m_b*IT_0254;
    const ccomplex_t IT_0256 = mtylt_C0iC(15, IT_0019, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0020, IT_0021, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0257 = IT_0027*IT_0256;
    const ccomplex_t IT_0258 = m_s*IT_0257;
    const ccomplex_t IT_0259 = m_b*IT_0030;
    const ccomplex_t IT_0260 = -IT_0252 + -IT_0255 + -IT_0258 + -IT_0259;
    const ccomplex_t IT_0261 = IT_0250 + IT_0260;
    const ccomplex_t IT_0262 = IT_0244*IT_0261;
    const ccomplex_t IT_0263 = IT_0037*IT_0076;
    const ccomplex_t IT_0264 = 0.101321183642338*IT_0263;
    const ccomplex_t IT_0265 = IT_0221*IT_0264;
    const ccomplex_t IT_0266 = IT_0061*IT_0083;
    const ccomplex_t IT_0267 = 0.101321183642338*IT_0266;
    const ccomplex_t IT_0268 = IT_0241*IT_0267;
    const ccomplex_t IT_0269 = IT_0016*IT_0069;
    const ccomplex_t IT_0270 = 0.101321183642338*IT_0269;
    const ccomplex_t IT_0271 = IT_0261*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0012;
    const ccomplex_t IT_0273 = 0.5*IT_0272;
    const ccomplex_t IT_0274 = IT_0189*IT_0273;
    const ccomplex_t IT_0275 = 0.101321183642338*IT_0274;
    const ccomplex_t IT_0276 = mtylt_C0iC(3, IT_0020, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0019, IT_0041, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0277 = 2*IT_0024;
    const ccomplex_t IT_0278 = IT_0276*IT_0277;
    const ccomplex_t IT_0279 = m_s*IT_0278;
    const ccomplex_t IT_0280 = IT_0097*IT_0277;
    const ccomplex_t IT_0281 = m_s*IT_0280;
    const ccomplex_t IT_0282 = IT_0025*IT_0097;
    const ccomplex_t IT_0283 = m_b*IT_0282;
    const ccomplex_t IT_0284 = IT_0025*IT_0276;
    const ccomplex_t IT_0285 = m_s*IT_0284;
    const ccomplex_t IT_0286 = mtylt_C0iC(15, IT_0020, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0019, IT_0041, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0287 = IT_0027*IT_0286;
    const ccomplex_t IT_0288 = m_s*IT_0287;
    const ccomplex_t IT_0289 = mtylt_C0iC(18, IT_0020, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0019, IT_0041, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0290 = IT_0027*IT_0289;
    const ccomplex_t IT_0291 = m_b*IT_0290;
    const ccomplex_t IT_0292 = IT_0279 + IT_0281 + IT_0283 + IT_0285 + IT_0288
       + IT_0291;
    const ccomplex_t IT_0293 = IT_0024*IT_0097;
    const ccomplex_t IT_0294 = m_s*IT_0293;
    const ccomplex_t IT_0295 = -IT_0294;
    const ccomplex_t IT_0296 = IT_0292 + IT_0295;
    const ccomplex_t IT_0297 = IT_0275*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0012;
    const ccomplex_t IT_0299 = 0.5*IT_0298;
    const ccomplex_t IT_0300 = IT_0194*IT_0299;
    const ccomplex_t IT_0301 = 0.101321183642338*IT_0300;
    const ccomplex_t IT_0302 = mtylt_C0iC(3, IT_0020, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0019, IT_0050, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0303 = IT_0025*IT_0302;
    const ccomplex_t IT_0304 = m_s*IT_0303;
    const ccomplex_t IT_0305 = IT_0025*IT_0111;
    const ccomplex_t IT_0306 = m_b*IT_0305;
    const ccomplex_t IT_0307 = IT_0277*IT_0302;
    const ccomplex_t IT_0308 = m_s*IT_0307;
    const ccomplex_t IT_0309 = IT_0111*IT_0277;
    const ccomplex_t IT_0310 = m_s*IT_0309;
    const ccomplex_t IT_0311 = mtylt_C0iC(15, IT_0020, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0019, IT_0050, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0312 = IT_0027*IT_0311;
    const ccomplex_t IT_0313 = m_s*IT_0312;
    const ccomplex_t IT_0314 = mtylt_C0iC(18, IT_0020, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0019, IT_0050, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0315 = IT_0027*IT_0314;
    const ccomplex_t IT_0316 = m_b*IT_0315;
    const ccomplex_t IT_0317 = IT_0304 + IT_0306 + IT_0308 + IT_0310 + IT_0313
       + IT_0316;
    const ccomplex_t IT_0318 = IT_0024*IT_0111;
    const ccomplex_t IT_0319 = m_s*IT_0318;
    const ccomplex_t IT_0320 = -IT_0319;
    const ccomplex_t IT_0321 = IT_0317 + IT_0320;
    const ccomplex_t IT_0322 = IT_0301*IT_0321;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0012;
    const ccomplex_t IT_0324 = 0.5*IT_0323;
    const ccomplex_t IT_0325 = IT_0199*IT_0324;
    const ccomplex_t IT_0326 = 0.101321183642338*IT_0325;
    const ccomplex_t IT_0327 = mtylt_C0iC(3, IT_0020, IT_0019 + IT_0020 + (-2)
      *s_12, IT_0019, IT_0021, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0328 = IT_0025*IT_0327;
    const ccomplex_t IT_0329 = m_s*IT_0328;
    const ccomplex_t IT_0330 = IT_0277*IT_0327;
    const ccomplex_t IT_0331 = m_s*IT_0330;
    const ccomplex_t IT_0332 = IT_0125*IT_0277;
    const ccomplex_t IT_0333 = m_s*IT_0332;
    const ccomplex_t IT_0334 = mtylt_C0iC(15, IT_0020, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0019, IT_0021, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0335 = IT_0027*IT_0334;
    const ccomplex_t IT_0336 = m_s*IT_0335;
    const ccomplex_t IT_0337 = IT_0025*IT_0125;
    const ccomplex_t IT_0338 = m_b*IT_0337;
    const ccomplex_t IT_0339 = mtylt_C0iC(18, IT_0020, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0019, IT_0021, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0340 = IT_0027*IT_0339;
    const ccomplex_t IT_0341 = m_b*IT_0340;
    const ccomplex_t IT_0342 = IT_0329 + IT_0331 + IT_0333 + IT_0336 + IT_0338
       + IT_0341;
    const ccomplex_t IT_0343 = IT_0024*IT_0125;
    const ccomplex_t IT_0344 = m_s*IT_0343;
    const ccomplex_t IT_0345 = -IT_0344;
    const ccomplex_t IT_0346 = IT_0342 + IT_0345;
    const ccomplex_t IT_0347 = IT_0326*IT_0346;
    const ccomplex_t IT_0348 = IT_0156*IT_0189*IT_0273;
    const ccomplex_t IT_0349 = 0.101321183642338*IT_0348;
    const ccomplex_t IT_0350 = m_s*IT_0161;
    const ccomplex_t IT_0351 = mtylt_C0iC(15, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0041, IT_0093, IT_0041, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0352 = m_s*IT_0351;
    const ccomplex_t IT_0353 = mtylt_C0iC(18, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0041, IT_0093, IT_0041, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0354 = m_s*IT_0353;
    const ccomplex_t IT_0355 = m_s*IT_0160;
    const ccomplex_t IT_0356 = IT_0350 + IT_0352 + IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = m_b*IT_0161;
    const ccomplex_t IT_0358 = mtylt_C0iC(12, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0041, IT_0093, IT_0041, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0359 = m_b*IT_0358;
    const ccomplex_t IT_0360 = m_b*IT_0351;
    const ccomplex_t IT_0361 = m_b*IT_0353;
    const ccomplex_t IT_0362 = m_b*IT_0160;
    const ccomplex_t IT_0363 = -IT_0357 + -IT_0359 + (-2)*IT_0360 + -IT_0361 +
       -IT_0362;
    const ccomplex_t IT_0364 = IT_0356 + IT_0363;
    const ccomplex_t IT_0365 = IT_0349*IT_0364;
    const ccomplex_t IT_0366 = IT_0035*IT_0037*IT_0156;
    const ccomplex_t IT_0367 = IT_0039*IT_0366;
    const ccomplex_t IT_0368 = IT_0155*IT_0367;
    const ccomplex_t IT_0369 = IT_0035*IT_0074*IT_0156;
    const ccomplex_t IT_0370 = 0.101321183642338*IT_0369;
    const ccomplex_t IT_0371 = m_s*IT_0153;
    const ccomplex_t IT_0372 = mtylt_C0iC(15, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0041, IT_0041, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0373 = m_s*IT_0372;
    const ccomplex_t IT_0374 = mtylt_C0iC(12, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0041, IT_0041, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0375 = m_s*IT_0374;
    const ccomplex_t IT_0376 = IT_0371 + IT_0373 + IT_0375;
    const ccomplex_t IT_0377 = m_b*IT_0153;
    const ccomplex_t IT_0378 = m_b*IT_0154;
    const ccomplex_t IT_0379 = mtylt_C0iC(18, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0041, IT_0041, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0380 = m_b*IT_0379;
    const ccomplex_t IT_0381 = m_b*IT_0374;
    const ccomplex_t IT_0382 = m_b*IT_0372;
    const ccomplex_t IT_0383 = -IT_0377 + -IT_0378 + -IT_0380 + -IT_0381 + (-2
      )*IT_0382;
    const ccomplex_t IT_0384 = IT_0376 + IT_0383;
    const ccomplex_t IT_0385 = IT_0370*IT_0384;
    const ccomplex_t IT_0386 = IT_0088*IT_0090*IT_0156;
    const ccomplex_t IT_0387 = IT_0039*IT_0386;
    const ccomplex_t IT_0388 = IT_0163*IT_0387;
    const ccomplex_t IT_0389 = IT_0088*IT_0131*IT_0156;
    const ccomplex_t IT_0390 = 0.101321183642338*IT_0389;
    const ccomplex_t IT_0391 = IT_0350 + IT_0352 + IT_0354;
    const ccomplex_t IT_0392 = IT_0363 + IT_0391;
    const ccomplex_t IT_0393 = IT_0390*IT_0392;
    const ccomplex_t IT_0394 = IT_0037*IT_0076*IT_0156;
    const ccomplex_t IT_0395 = 0.101321183642338*IT_0394;
    const ccomplex_t IT_0396 = IT_0384*IT_0395;
    const ccomplex_t IT_0397 = IT_0090*IT_0133*IT_0156;
    const ccomplex_t IT_0398 = 0.101321183642338*IT_0397;
    const ccomplex_t IT_0399 = IT_0392*IT_0398;
    const ccomplex_t IT_0400 = IT_0156*IT_0194*IT_0299;
    const ccomplex_t IT_0401 = 0.101321183642338*IT_0400;
    const ccomplex_t IT_0402 = mtylt_C0iC(18, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0050, IT_0093, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0403 = m_s*IT_0402;
    const ccomplex_t IT_0404 = m_s*IT_0174;
    const ccomplex_t IT_0405 = m_s*IT_0175;
    const ccomplex_t IT_0406 = mtylt_C0iC(15, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0050, IT_0093, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0407 = m_s*IT_0406;
    const ccomplex_t IT_0408 = IT_0403 + IT_0404 + IT_0405 + IT_0407;
    const ccomplex_t IT_0409 = mtylt_C0iC(12, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0050, IT_0093, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0410 = m_b*IT_0409;
    const ccomplex_t IT_0411 = m_b*IT_0174;
    const ccomplex_t IT_0412 = m_b*IT_0175;
    const ccomplex_t IT_0413 = m_b*IT_0406;
    const ccomplex_t IT_0414 = m_b*IT_0402;
    const ccomplex_t IT_0415 = -IT_0410 + -IT_0411 + -IT_0412 + (-2)*IT_0413 +
       -IT_0414;
    const ccomplex_t IT_0416 = IT_0408 + IT_0415;
    const ccomplex_t IT_0417 = IT_0401*IT_0416;
    const ccomplex_t IT_0418 = IT_0059*IT_0081*IT_0156;
    const ccomplex_t IT_0419 = 0.101321183642338*IT_0418;
    const ccomplex_t IT_0420 = m_s*IT_0170;
    const ccomplex_t IT_0421 = mtylt_C0iC(15, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0050, IT_0050, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0422 = m_s*IT_0421;
    const ccomplex_t IT_0423 = mtylt_C0iC(12, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0050, IT_0050, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0424 = m_s*IT_0423;
    const ccomplex_t IT_0425 = IT_0420 + IT_0422 + IT_0424;
    const ccomplex_t IT_0426 = m_b*IT_0170;
    const ccomplex_t IT_0427 = m_b*IT_0171;
    const ccomplex_t IT_0428 = m_b*IT_0423;
    const ccomplex_t IT_0429 = m_b*IT_0421;
    const ccomplex_t IT_0430 = mtylt_C0iC(18, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0050, IT_0050, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0431 = m_b*IT_0430;
    const ccomplex_t IT_0432 = -IT_0426 + -IT_0427 + -IT_0428 + (-2)*IT_0429 +
       -IT_0431;
    const ccomplex_t IT_0433 = IT_0425 + IT_0432;
    const ccomplex_t IT_0434 = IT_0419*IT_0433;
    const ccomplex_t IT_0435 = IT_0103*IT_0138*IT_0156;
    const ccomplex_t IT_0436 = 0.101321183642338*IT_0435;
    const ccomplex_t IT_0437 = IT_0403 + IT_0405 + IT_0407;
    const ccomplex_t IT_0438 = IT_0415 + IT_0437;
    const ccomplex_t IT_0439 = IT_0436*IT_0438;
    const ccomplex_t IT_0440 = IT_0061*IT_0083*IT_0156;
    const ccomplex_t IT_0441 = 0.101321183642338*IT_0440;
    const ccomplex_t IT_0442 = IT_0433*IT_0441;
    const ccomplex_t IT_0443 = IT_0105*IT_0140*IT_0156;
    const ccomplex_t IT_0444 = 0.101321183642338*IT_0443;
    const ccomplex_t IT_0445 = IT_0438*IT_0444;
    const ccomplex_t IT_0446 = IT_0156*IT_0199*IT_0324;
    const ccomplex_t IT_0447 = 0.101321183642338*IT_0446;
    const ccomplex_t IT_0448 = mtylt_C0iC(6, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0021, IT_0093, IT_0021, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0449 = m_s*IT_0448;
    const ccomplex_t IT_0450 = mtylt_C0iC(18, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0021, IT_0093, IT_0021, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0451 = m_s*IT_0450;
    const ccomplex_t IT_0452 = mtylt_C0iC(15, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0021, IT_0093, IT_0021, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0453 = m_s*IT_0452;
    const ccomplex_t IT_0454 = mtylt_C0iC(3, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0021, IT_0093, IT_0021, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0455 = m_s*IT_0454;
    const ccomplex_t IT_0456 = IT_0449 + IT_0451 + IT_0453 + IT_0455;
    const ccomplex_t IT_0457 = m_b*IT_0448;
    const ccomplex_t IT_0458 = m_b*IT_0452;
    const ccomplex_t IT_0459 = m_b*IT_0450;
    const ccomplex_t IT_0460 = mtylt_C0iC(12, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0021, IT_0093, IT_0021, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0461 = m_b*IT_0460;
    const ccomplex_t IT_0462 = m_b*IT_0454;
    const ccomplex_t IT_0463 = -IT_0457 + (-2)*IT_0458 + -IT_0459 + -IT_0461 +
       -IT_0462;
    const ccomplex_t IT_0464 = IT_0456 + IT_0463;
    const ccomplex_t IT_0465 = IT_0447*IT_0464;
    const ccomplex_t IT_0466 = IT_0067*IT_0069*IT_0156;
    const ccomplex_t IT_0467 = IT_0006*IT_0466;
    const ccomplex_t IT_0468 = mtylt_C0iC(0, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0021, IT_0021, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0469 = mtylt_C0iC(3, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0021, IT_0021, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0470 = mtylt_C0iC(6, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0021, IT_0021, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0471 = IT_0468 + IT_0469 + IT_0470;
    const ccomplex_t IT_0472 = IT_0467*IT_0471;
    const ccomplex_t IT_0473 = IT_0014*IT_0067*IT_0156;
    const ccomplex_t IT_0474 = 0.101321183642338*IT_0473;
    const ccomplex_t IT_0475 = m_s*IT_0469;
    const ccomplex_t IT_0476 = mtylt_C0iC(12, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0021, IT_0021, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0477 = m_s*IT_0476;
    const ccomplex_t IT_0478 = mtylt_C0iC(15, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0021, IT_0021, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0479 = m_s*IT_0478;
    const ccomplex_t IT_0480 = IT_0475 + IT_0477 + IT_0479;
    const ccomplex_t IT_0481 = m_b*IT_0469;
    const ccomplex_t IT_0482 = m_b*IT_0476;
    const ccomplex_t IT_0483 = m_b*IT_0478;
    const ccomplex_t IT_0484 = m_b*IT_0470;
    const ccomplex_t IT_0485 = mtylt_C0iC(18, IT_0019 + IT_0020 + (-2)*s_12,
       IT_0020, IT_0019, IT_0021, IT_0021, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0486 = m_b*IT_0485;
    const ccomplex_t IT_0487 = -IT_0481 + -IT_0482 + (-2)*IT_0483 + -IT_0484 +
       -IT_0486;
    const ccomplex_t IT_0488 = IT_0480 + IT_0487;
    const ccomplex_t IT_0489 = IT_0474*IT_0488;
    const ccomplex_t IT_0490 = IT_0117*IT_0119*IT_0156;
    const ccomplex_t IT_0491 = IT_0006*IT_0490;
    const ccomplex_t IT_0492 = mtylt_C0iC(0, IT_0019, IT_0020, IT_0019 +
       IT_0020 + (-2)*s_12, IT_0021, IT_0093, IT_0021, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0493 = IT_0448 + IT_0454 + IT_0492;
    const ccomplex_t IT_0494 = IT_0491*IT_0493;
    const ccomplex_t IT_0495 = IT_0117*IT_0145*IT_0156;
    const ccomplex_t IT_0496 = 0.101321183642338*IT_0495;
    const ccomplex_t IT_0497 = IT_0449 + IT_0451 + IT_0453;
    const ccomplex_t IT_0498 = IT_0463 + IT_0497;
    const ccomplex_t IT_0499 = IT_0496*IT_0498;
    const ccomplex_t IT_0500 = IT_0016*IT_0069*IT_0156;
    const ccomplex_t IT_0501 = 0.101321183642338*IT_0500;
    const ccomplex_t IT_0502 = IT_0488*IT_0501;
    const ccomplex_t IT_0503 = IT_0014*IT_0016*IT_0156;
    const ccomplex_t IT_0504 = IT_0006*IT_0503;
    const ccomplex_t IT_0505 = IT_0471*IT_0504;
    const ccomplex_t IT_0506 = IT_0119*IT_0147*IT_0156;
    const ccomplex_t IT_0507 = 0.101321183642338*IT_0506;
    const ccomplex_t IT_0508 = IT_0498*IT_0507;
    const ccomplex_t IT_0509 = IT_0145*IT_0147*IT_0156;
    const ccomplex_t IT_0510 = IT_0006*IT_0509;
    const ccomplex_t IT_0511 = IT_0493*IT_0510;
    const ccomplex_t IT_0512 = IT_0088*IT_0131;
    const ccomplex_t IT_0513 = 0.101321183642338*IT_0512;
    const ccomplex_t IT_0514 = IT_0283 + IT_0285;
    const ccomplex_t IT_0515 = m_b*IT_0098;
    const ccomplex_t IT_0516 = IT_0027*IT_0276;
    const ccomplex_t IT_0517 = m_s*IT_0516;
    const ccomplex_t IT_0518 = -IT_0288 + -IT_0291 + -IT_0515 + -IT_0517;
    const ccomplex_t IT_0519 = IT_0514 + IT_0518;
    const ccomplex_t IT_0520 = IT_0513*IT_0519;
    const ccomplex_t IT_0521 = IT_0103*IT_0138;
    const ccomplex_t IT_0522 = 0.101321183642338*IT_0521;
    const ccomplex_t IT_0523 = IT_0304 + IT_0306;
    const ccomplex_t IT_0524 = m_b*IT_0112;
    const ccomplex_t IT_0525 = IT_0027*IT_0302;
    const ccomplex_t IT_0526 = m_s*IT_0525;
    const ccomplex_t IT_0527 = -IT_0313 + -IT_0316 + -IT_0524 + -IT_0526;
    const ccomplex_t IT_0528 = IT_0523 + IT_0527;
    const ccomplex_t IT_0529 = IT_0522*IT_0528;
    const ccomplex_t IT_0530 = IT_0117*IT_0145;
    const ccomplex_t IT_0531 = 0.101321183642338*IT_0530;
    const ccomplex_t IT_0532 = IT_0329 + IT_0338;
    const ccomplex_t IT_0533 = IT_0027*IT_0327;
    const ccomplex_t IT_0534 = m_s*IT_0533;
    const ccomplex_t IT_0535 = m_b*IT_0126;
    const ccomplex_t IT_0536 = -IT_0336 + -IT_0341 + -IT_0534 + -IT_0535;
    const ccomplex_t IT_0537 = IT_0532 + IT_0536;
    const ccomplex_t IT_0538 = IT_0531*IT_0537;
    const ccomplex_t IT_0539 = IT_0090*IT_0133;
    const ccomplex_t IT_0540 = 0.101321183642338*IT_0539;
    const ccomplex_t IT_0541 = IT_0519*IT_0540;
    const ccomplex_t IT_0542 = IT_0105*IT_0140;
    const ccomplex_t IT_0543 = 0.101321183642338*IT_0542;
    const ccomplex_t IT_0544 = IT_0528*IT_0543;
    const ccomplex_t IT_0545 = IT_0119*IT_0147;
    const ccomplex_t IT_0546 = 0.101321183642338*IT_0545;
    const ccomplex_t IT_0547 = IT_0537*IT_0546;
    const ccomplex_t IT_0548 = 2*IT_0159 + 2*IT_0166 + 2*IT_0173 + 2*IT_0180 +
       2*IT_0183 + 2*IT_0186 + (-2)*IT_0192 + (-2)*IT_0197 + (-2)*IT_0202 + 
      -IT_0222 + -IT_0242 + -IT_0262 + -IT_0265 + -IT_0268 + -IT_0271 + 2
      *IT_0297 + 2*IT_0322 + 2*IT_0347 + 4*IT_0365 + 2*IT_0368 + 2*IT_0385 + 2
      *IT_0388 + 2*IT_0393 + 2*IT_0396 + 2*IT_0399 + 4*IT_0417 + 2*IT_0434 + 2
      *IT_0439 + 2*IT_0442 + 2*IT_0445 + 4*IT_0465 + 2*IT_0472 + 2*IT_0489 + 2
      *IT_0494 + 2*IT_0499 + 2*IT_0502 + 2*IT_0505 + 2*IT_0508 + 2*IT_0511 + 
      -IT_0520 + -IT_0529 + -IT_0538 + -IT_0541 + -IT_0544 + -IT_0547;
    const ccomplex_t IT_0549 = IT_0151 + IT_0548;
    const ccomplex_t IT_0550 = IT_0005*IT_0549;
    const ccomplex_t IT_0551 = 1.4142135623731*IT_0550;
    const ccomplex_t IT_0552 = 0.125*IT_0551;
    const ccomplex_t IT_0553 = IT_0109 + IT_0525;
    const ccomplex_t IT_0554 = IT_0107*IT_0553;
    const ccomplex_t IT_0555 = IT_0123 + IT_0533;
    const ccomplex_t IT_0556 = IT_0121*IT_0555;
    const ccomplex_t IT_0557 = IT_0043 + IT_0211;
    const ccomplex_t IT_0558 = IT_0040*IT_0557;
    const ccomplex_t IT_0559 = IT_0052 + IT_0231;
    const ccomplex_t IT_0560 = IT_0064*IT_0559;
    const ccomplex_t IT_0561 = IT_0026 + IT_0251;
    const ccomplex_t IT_0562 = IT_0071*IT_0561;
    const ccomplex_t IT_0563 = IT_0078*IT_0557;
    const ccomplex_t IT_0564 = IT_0085*IT_0559;
    const ccomplex_t IT_0565 = IT_0018*IT_0561;
    const ccomplex_t IT_0566 = IT_0095 + IT_0516;
    const ccomplex_t IT_0567 = IT_0092*IT_0566;
    const ccomplex_t IT_0568 = IT_0135*IT_0566;
    const ccomplex_t IT_0569 = IT_0142*IT_0553;
    const ccomplex_t IT_0570 = IT_0149*IT_0555;
    const ccomplex_t IT_0571 = IT_0554 + IT_0556 + IT_0558 + IT_0560 + IT_0562
       + IT_0563 + IT_0564 + IT_0565 + IT_0567 + IT_0568 + IT_0569 + IT_0570;
    const ccomplex_t IT_0572 = IT_0168*IT_0170;
    const ccomplex_t IT_0573 = IT_0469*IT_0504;
    const ccomplex_t IT_0574 = IT_0448*IT_0510;
    const ccomplex_t IT_0575 = IT_0133*IT_0187*IT_0273;
    const ccomplex_t IT_0576 = 0.101321183642338*IT_0575;
    const ccomplex_t IT_0577 = IT_0276*IT_0576;
    const ccomplex_t IT_0578 = IT_0140*IT_0187*IT_0299;
    const ccomplex_t IT_0579 = 0.101321183642338*IT_0578;
    const ccomplex_t IT_0580 = IT_0302*IT_0579;
    const ccomplex_t IT_0581 = IT_0147*IT_0187*IT_0324;
    const ccomplex_t IT_0582 = 0.101321183642338*IT_0581;
    const ccomplex_t IT_0583 = IT_0327*IT_0582;
    const ccomplex_t IT_0584 = IT_0161*IT_0165;
    const ccomplex_t IT_0585 = IT_0170*IT_0182;
    const ccomplex_t IT_0586 = IT_0448*IT_0491;
    const ccomplex_t IT_0587 = mtylt_C0iC(18, IT_0019, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0020, IT_0041, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0588 = IT_0027*IT_0587;
    const ccomplex_t IT_0589 = m_s*IT_0588;
    const ccomplex_t IT_0590 = m_b*IT_0214;
    const ccomplex_t IT_0591 = IT_0207 + IT_0209 + IT_0589 + IT_0590;
    const ccomplex_t IT_0592 = IT_0204*IT_0591;
    const ccomplex_t IT_0593 = mtylt_C0iC(18, IT_0019, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0020, IT_0050, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0594 = IT_0027*IT_0593;
    const ccomplex_t IT_0595 = m_s*IT_0594;
    const ccomplex_t IT_0596 = m_b*IT_0234;
    const ccomplex_t IT_0597 = IT_0227 + IT_0229 + IT_0595 + IT_0596;
    const ccomplex_t IT_0598 = IT_0224*IT_0597;
    const ccomplex_t IT_0599 = mtylt_C0iC(18, IT_0019, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0020, IT_0021, IT_0022, IT_0022, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0600 = IT_0027*IT_0599;
    const ccomplex_t IT_0601 = m_s*IT_0600;
    const ccomplex_t IT_0602 = m_b*IT_0257;
    const ccomplex_t IT_0603 = IT_0247 + IT_0249 + IT_0601 + IT_0602;
    const ccomplex_t IT_0604 = IT_0244*IT_0603;
    const ccomplex_t IT_0605 = IT_0264*IT_0591;
    const ccomplex_t IT_0606 = IT_0267*IT_0597;
    const ccomplex_t IT_0607 = IT_0270*IT_0603;
    const ccomplex_t IT_0608 = m_b*IT_0284;
    const ccomplex_t IT_0609 = IT_0285 + IT_0608;
    const ccomplex_t IT_0610 = m_b*IT_0293;
    const ccomplex_t IT_0611 = mtylt_C0iC(12, IT_0020, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0019, IT_0041, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0612 = IT_0027*IT_0611;
    const ccomplex_t IT_0613 = m_s*IT_0612;
    const ccomplex_t IT_0614 = m_b*IT_0287;
    const ccomplex_t IT_0615 = 2*IT_0283 + -IT_0515 + -IT_0517 + -IT_0610 + 
      -IT_0613 + -IT_0614;
    const ccomplex_t IT_0616 = IT_0609 + IT_0615;
    const ccomplex_t IT_0617 = IT_0275*IT_0616;
    const ccomplex_t IT_0618 = m_b*IT_0303;
    const ccomplex_t IT_0619 = IT_0304 + IT_0618;
    const ccomplex_t IT_0620 = mtylt_C0iC(12, IT_0020, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0019, IT_0050, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0621 = IT_0027*IT_0620;
    const ccomplex_t IT_0622 = m_s*IT_0621;
    const ccomplex_t IT_0623 = m_b*IT_0318;
    const ccomplex_t IT_0624 = m_b*IT_0312;
    const ccomplex_t IT_0625 = 2*IT_0306 + -IT_0524 + -IT_0526 + -IT_0622 + 
      -IT_0623 + -IT_0624;
    const ccomplex_t IT_0626 = IT_0619 + IT_0625;
    const ccomplex_t IT_0627 = IT_0301*IT_0626;
    const ccomplex_t IT_0628 = m_b*IT_0328;
    const ccomplex_t IT_0629 = IT_0329 + IT_0628;
    const ccomplex_t IT_0630 = m_b*IT_0335;
    const ccomplex_t IT_0631 = mtylt_C0iC(12, IT_0020, IT_0019 + IT_0020 + (-2
      )*s_12, IT_0019, IT_0021, IT_0093, IT_0093, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0632 = IT_0027*IT_0631;
    const ccomplex_t IT_0633 = m_s*IT_0632;
    const ccomplex_t IT_0634 = m_b*IT_0343;
    const ccomplex_t IT_0635 = 2*IT_0338 + -IT_0534 + -IT_0535 + -IT_0630 + 
      -IT_0633 + -IT_0634;
    const ccomplex_t IT_0636 = IT_0629 + IT_0635;
    const ccomplex_t IT_0637 = IT_0326*IT_0636;
    const ccomplex_t IT_0638 = IT_0350 + IT_0354;
    const ccomplex_t IT_0639 = -IT_0357 + -IT_0360 + -IT_0361 + -IT_0362;
    const ccomplex_t IT_0640 = IT_0638 + IT_0639;
    const ccomplex_t IT_0641 = IT_0349*IT_0640;
    const ccomplex_t IT_0642 = IT_0153*IT_0367;
    const ccomplex_t IT_0643 = IT_0371 + IT_0375;
    const ccomplex_t IT_0644 = -IT_0377 + -IT_0381 + -IT_0382;
    const ccomplex_t IT_0645 = IT_0643 + IT_0644;
    const ccomplex_t IT_0646 = IT_0370*IT_0645;
    const ccomplex_t IT_0647 = IT_0161*IT_0387;
    const ccomplex_t IT_0648 = -IT_0357 + -IT_0360 + -IT_0361;
    const ccomplex_t IT_0649 = IT_0638 + IT_0648;
    const ccomplex_t IT_0650 = IT_0390*IT_0649;
    const ccomplex_t IT_0651 = IT_0395*IT_0645;
    const ccomplex_t IT_0652 = IT_0153*IT_0158;
    const ccomplex_t IT_0653 = IT_0398*IT_0649;
    const ccomplex_t IT_0654 = IT_0403 + IT_0405;
    const ccomplex_t IT_0655 = -IT_0411 + -IT_0412 + -IT_0413 + -IT_0414;
    const ccomplex_t IT_0656 = IT_0654 + IT_0655;
    const ccomplex_t IT_0657 = IT_0401*IT_0656;
    const ccomplex_t IT_0658 = IT_0420 + IT_0424;
    const ccomplex_t IT_0659 = -IT_0426 + -IT_0428 + -IT_0429;
    const ccomplex_t IT_0660 = IT_0658 + IT_0659;
    const ccomplex_t IT_0661 = IT_0419*IT_0660;
    const ccomplex_t IT_0662 = IT_0175*IT_0179;
    const ccomplex_t IT_0663 = -IT_0412 + -IT_0413 + -IT_0414;
    const ccomplex_t IT_0664 = IT_0654 + IT_0663;
    const ccomplex_t IT_0665 = IT_0436*IT_0664;
    const ccomplex_t IT_0666 = IT_0441*IT_0660;
    const ccomplex_t IT_0667 = IT_0444*IT_0664;
    const ccomplex_t IT_0668 = IT_0175*IT_0185;
    const ccomplex_t IT_0669 = IT_0449 + IT_0451;
    const ccomplex_t IT_0670 = -IT_0457 + -IT_0458 + -IT_0459 + -IT_0462;
    const ccomplex_t IT_0671 = IT_0669 + IT_0670;
    const ccomplex_t IT_0672 = IT_0447*IT_0671;
    const ccomplex_t IT_0673 = IT_0467*IT_0469;
    const ccomplex_t IT_0674 = IT_0475 + IT_0477;
    const ccomplex_t IT_0675 = -IT_0481 + -IT_0482 + -IT_0483;
    const ccomplex_t IT_0676 = IT_0674 + IT_0675;
    const ccomplex_t IT_0677 = IT_0474*IT_0676;
    const ccomplex_t IT_0678 = -IT_0457 + -IT_0458 + -IT_0459;
    const ccomplex_t IT_0679 = IT_0669 + IT_0678;
    const ccomplex_t IT_0680 = IT_0496*IT_0679;
    const ccomplex_t IT_0681 = IT_0501*IT_0676;
    const ccomplex_t IT_0682 = IT_0507*IT_0679;
    const ccomplex_t IT_0683 = IT_0283 + IT_0285 + IT_0613 + IT_0614;
    const ccomplex_t IT_0684 = IT_0513*IT_0683;
    const ccomplex_t IT_0685 = IT_0304 + IT_0306 + IT_0622 + IT_0624;
    const ccomplex_t IT_0686 = IT_0522*IT_0685;
    const ccomplex_t IT_0687 = IT_0329 + IT_0338 + IT_0630 + IT_0633;
    const ccomplex_t IT_0688 = IT_0531*IT_0687;
    const ccomplex_t IT_0689 = IT_0540*IT_0683;
    const ccomplex_t IT_0690 = IT_0543*IT_0685;
    const ccomplex_t IT_0691 = IT_0546*IT_0687;
    const ccomplex_t IT_0692 = 2*IT_0572 + 2*IT_0573 + 2*IT_0574 + 2*IT_0577 +
       2*IT_0580 + 2*IT_0583 + 2*IT_0584 + 2*IT_0585 + 2*IT_0586 + -IT_0592 + 
      -IT_0598 + -IT_0604 + -IT_0605 + -IT_0606 + -IT_0607 + 2*IT_0617 + 2
      *IT_0627 + 2*IT_0637 + 4*IT_0641 + 2*IT_0642 + 2*IT_0646 + 2*IT_0647 + 2
      *IT_0650 + 2*IT_0651 + 2*IT_0652 + 2*IT_0653 + 4*IT_0657 + 2*IT_0661 + 2
      *IT_0662 + 2*IT_0665 + 2*IT_0666 + 2*IT_0667 + 2*IT_0668 + 4*IT_0672 + 2
      *IT_0673 + 2*IT_0677 + 2*IT_0680 + 2*IT_0681 + 2*IT_0682 + -IT_0684 + 
      -IT_0686 + -IT_0688 + -IT_0689 + -IT_0690 + -IT_0691;
    const ccomplex_t IT_0693 = IT_0571 + IT_0692;
    const ccomplex_t IT_0694 = IT_0005*IT_0693;
    const ccomplex_t IT_0695 = 1.4142135623731*IT_0694;
    const ccomplex_t IT_0696 = (-0.125)*IT_0695;
    const ccomplex_t IT_0697 = IT_0207 + IT_0216 + IT_0219;
    const ccomplex_t IT_0698 = -IT_0209 + -IT_0212 + -IT_0215;
    const ccomplex_t IT_0699 = IT_0697 + IT_0698;
    const ccomplex_t IT_0700 = IT_0264*IT_0699;
    const ccomplex_t IT_0701 = IT_0227 + IT_0236 + IT_0239;
    const ccomplex_t IT_0702 = -IT_0229 + -IT_0232 + -IT_0235;
    const ccomplex_t IT_0703 = IT_0701 + IT_0702;
    const ccomplex_t IT_0704 = IT_0267*IT_0703;
    const ccomplex_t IT_0705 = -IT_0249 + -IT_0252 + -IT_0258;
    const ccomplex_t IT_0706 = IT_0247 + IT_0255 + IT_0259;
    const ccomplex_t IT_0707 = IT_0705 + IT_0706;
    const ccomplex_t IT_0708 = IT_0270*IT_0707;
    const ccomplex_t IT_0709 = IT_0285 + IT_0291 + IT_0515;
    const ccomplex_t IT_0710 = -IT_0283 + -IT_0288 + -IT_0517;
    const ccomplex_t IT_0711 = IT_0709 + IT_0710;
    const ccomplex_t IT_0712 = IT_0540*IT_0711;
    const ccomplex_t IT_0713 = IT_0304 + IT_0316 + IT_0524;
    const ccomplex_t IT_0714 = -IT_0306 + -IT_0313 + -IT_0526;
    const ccomplex_t IT_0715 = IT_0713 + IT_0714;
    const ccomplex_t IT_0716 = IT_0543*IT_0715;
    const ccomplex_t IT_0717 = IT_0329 + IT_0341 + IT_0535;
    const ccomplex_t IT_0718 = -IT_0336 + -IT_0338 + -IT_0534;
    const ccomplex_t IT_0719 = IT_0717 + IT_0718;
    const ccomplex_t IT_0720 = IT_0546*IT_0719;
    const ccomplex_t IT_0721 = IT_0049 + IT_0065 + IT_0072 + IT_0101 + IT_0115
       + IT_0129 + IT_0700 + IT_0704 + IT_0708 + IT_0712 + IT_0716 + IT_0720;
    const ccomplex_t IT_0722 = IT_0204*IT_0699;
    const ccomplex_t IT_0723 = IT_0224*IT_0703;
    const ccomplex_t IT_0724 = IT_0244*IT_0707;
    const ccomplex_t IT_0725 = IT_0279 + IT_0281 + IT_0285 + IT_0288;
    const ccomplex_t IT_0726 = -IT_0283 + -IT_0291 + -IT_0294;
    const ccomplex_t IT_0727 = IT_0725 + IT_0726;
    const ccomplex_t IT_0728 = IT_0275*IT_0727;
    const ccomplex_t IT_0729 = IT_0304 + IT_0308 + IT_0310 + IT_0313;
    const ccomplex_t IT_0730 = -IT_0306 + -IT_0316 + -IT_0319;
    const ccomplex_t IT_0731 = IT_0729 + IT_0730;
    const ccomplex_t IT_0732 = IT_0301*IT_0731;
    const ccomplex_t IT_0733 = IT_0329 + IT_0331 + IT_0333 + IT_0336;
    const ccomplex_t IT_0734 = -IT_0338 + -IT_0341 + -IT_0344;
    const ccomplex_t IT_0735 = IT_0733 + IT_0734;
    const ccomplex_t IT_0736 = IT_0326*IT_0735;
    const ccomplex_t IT_0737 = 2*IT_0360;
    const ccomplex_t IT_0738 = IT_0350 + IT_0352 + IT_0354 + IT_0355 + IT_0357
       + IT_0359 + IT_0361 + IT_0362;
    const ccomplex_t IT_0739 = IT_0737 + IT_0738;
    const ccomplex_t IT_0740 = IT_0349*IT_0739;
    const ccomplex_t IT_0741 = IT_0371 + IT_0373 + IT_0375 + IT_0377 + IT_0378
       + IT_0380 + IT_0381;
    const ccomplex_t IT_0742 = 2*IT_0382;
    const ccomplex_t IT_0743 = IT_0741 + IT_0742;
    const ccomplex_t IT_0744 = IT_0370*IT_0743;
    const ccomplex_t IT_0745 = IT_0350 + IT_0352 + IT_0354 + IT_0357 + IT_0359
       + IT_0361 + IT_0362;
    const ccomplex_t IT_0746 = IT_0737 + IT_0745;
    const ccomplex_t IT_0747 = IT_0390*IT_0746;
    const ccomplex_t IT_0748 = IT_0395*IT_0743;
    const ccomplex_t IT_0749 = IT_0398*IT_0746;
    const ccomplex_t IT_0750 = 2*IT_0413;
    const ccomplex_t IT_0751 = IT_0403 + IT_0404 + IT_0405 + IT_0407 + IT_0410
       + IT_0411 + IT_0412 + IT_0414;
    const ccomplex_t IT_0752 = IT_0750 + IT_0751;
    const ccomplex_t IT_0753 = IT_0401*IT_0752;
    const ccomplex_t IT_0754 = IT_0420 + IT_0422 + IT_0424 + IT_0426 + IT_0427
       + IT_0428 + IT_0431;
    const ccomplex_t IT_0755 = 2*IT_0429;
    const ccomplex_t IT_0756 = IT_0754 + IT_0755;
    const ccomplex_t IT_0757 = IT_0419*IT_0756;
    const ccomplex_t IT_0758 = IT_0403 + IT_0405 + IT_0407 + IT_0410 + IT_0411
       + IT_0412 + IT_0414;
    const ccomplex_t IT_0759 = IT_0750 + IT_0758;
    const ccomplex_t IT_0760 = IT_0436*IT_0759;
    const ccomplex_t IT_0761 = IT_0441*IT_0756;
    const ccomplex_t IT_0762 = IT_0444*IT_0759;
    const ccomplex_t IT_0763 = IT_0449 + IT_0451 + IT_0453 + IT_0455 + IT_0457
       + IT_0459 + IT_0461 + IT_0462;
    const ccomplex_t IT_0764 = 2*IT_0458;
    const ccomplex_t IT_0765 = IT_0763 + IT_0764;
    const ccomplex_t IT_0766 = IT_0447*IT_0765;
    const ccomplex_t IT_0767 = IT_0475 + IT_0477 + IT_0479 + IT_0481 + IT_0482
       + IT_0484 + IT_0486;
    const ccomplex_t IT_0768 = 2*IT_0483;
    const ccomplex_t IT_0769 = IT_0767 + IT_0768;
    const ccomplex_t IT_0770 = IT_0474*IT_0769;
    const ccomplex_t IT_0771 = IT_0449 + IT_0451 + IT_0453 + IT_0457 + IT_0459
       + IT_0461 + IT_0462;
    const ccomplex_t IT_0772 = IT_0764 + IT_0771;
    const ccomplex_t IT_0773 = IT_0496*IT_0772;
    const ccomplex_t IT_0774 = IT_0501*IT_0769;
    const ccomplex_t IT_0775 = IT_0507*IT_0772;
    const ccomplex_t IT_0776 = IT_0513*IT_0711;
    const ccomplex_t IT_0777 = IT_0522*IT_0715;
    const ccomplex_t IT_0778 = IT_0531*IT_0719;
    const ccomplex_t IT_0779 = -IT_0033 + -IT_0079 + -IT_0086 + -IT_0136 + 
      -IT_0143 + -IT_0150 + (-2)*IT_0159 + (-2)*IT_0166 + 2*IT_0173 + 2*IT_0180 
      + (-2)*IT_0183 + (-2)*IT_0186 + (-2)*IT_0192 + (-2)*IT_0197 + (-2)*IT_0202
       + 2*IT_0368 + 2*IT_0388 + 2*IT_0472 + 2*IT_0494 + (-2)*IT_0505 + (-2)
      *IT_0511 + -IT_0722 + -IT_0723 + -IT_0724 + 2*IT_0728 + 2*IT_0732 + 2
      *IT_0736 + 4*IT_0740 + 2*IT_0744 + 2*IT_0747 + (-2)*IT_0748 + (-2)*IT_0749
       + 4*IT_0753 + 2*IT_0757 + 2*IT_0760 + (-2)*IT_0761 + (-2)*IT_0762 + 4
      *IT_0766 + 2*IT_0770 + 2*IT_0773 + (-2)*IT_0774 + (-2)*IT_0775 + -IT_0776 
      + -IT_0777 + -IT_0778;
    const ccomplex_t IT_0780 = IT_0721 + IT_0779;
    const ccomplex_t IT_0781 = IT_0005*IT_0780;
    const ccomplex_t IT_0782 = 1.4142135623731*IT_0781;
    const ccomplex_t IT_0783 = (-0.125)*IT_0782;
    const ccomplex_t IT_0784 = IT_0207 + IT_0589;
    const ccomplex_t IT_0785 = -IT_0209 + -IT_0590;
    const ccomplex_t IT_0786 = IT_0784 + IT_0785;
    const ccomplex_t IT_0787 = IT_0264*IT_0786;
    const ccomplex_t IT_0788 = -IT_0229 + -IT_0596;
    const ccomplex_t IT_0789 = IT_0227 + IT_0595;
    const ccomplex_t IT_0790 = IT_0788 + IT_0789;
    const ccomplex_t IT_0791 = IT_0267*IT_0790;
    const ccomplex_t IT_0792 = IT_0247 + IT_0601;
    const ccomplex_t IT_0793 = -IT_0249 + -IT_0602;
    const ccomplex_t IT_0794 = IT_0792 + IT_0793;
    const ccomplex_t IT_0795 = IT_0270*IT_0794;
    const ccomplex_t IT_0796 = -IT_0283 + -IT_0614;
    const ccomplex_t IT_0797 = IT_0285 + IT_0613;
    const ccomplex_t IT_0798 = IT_0796 + IT_0797;
    const ccomplex_t IT_0799 = IT_0540*IT_0798;
    const ccomplex_t IT_0800 = IT_0304 + IT_0622;
    const ccomplex_t IT_0801 = -IT_0306 + -IT_0624;
    const ccomplex_t IT_0802 = IT_0800 + IT_0801;
    const ccomplex_t IT_0803 = IT_0543*IT_0802;
    const ccomplex_t IT_0804 = IT_0329 + IT_0633;
    const ccomplex_t IT_0805 = -IT_0338 + -IT_0630;
    const ccomplex_t IT_0806 = IT_0804 + IT_0805;
    const ccomplex_t IT_0807 = IT_0546*IT_0806;
    const ccomplex_t IT_0808 = IT_0554 + IT_0556 + IT_0558 + IT_0560 + IT_0562
       + IT_0567 + IT_0787 + IT_0791 + IT_0795 + IT_0799 + IT_0803 + IT_0807;
    const ccomplex_t IT_0809 = IT_0204*IT_0786;
    const ccomplex_t IT_0810 = IT_0224*IT_0790;
    const ccomplex_t IT_0811 = IT_0244*IT_0794;
    const ccomplex_t IT_0812 = IT_0285 + IT_0515 + IT_0610 + IT_0614;
    const ccomplex_t IT_0813 = (-2)*IT_0283 + -IT_0517 + -IT_0608 + -IT_0613;
    const ccomplex_t IT_0814 = IT_0812 + IT_0813;
    const ccomplex_t IT_0815 = IT_0275*IT_0814;
    const ccomplex_t IT_0816 = IT_0304 + IT_0524 + IT_0623 + IT_0624;
    const ccomplex_t IT_0817 = (-2)*IT_0306 + -IT_0526 + -IT_0618 + -IT_0622;
    const ccomplex_t IT_0818 = IT_0816 + IT_0817;
    const ccomplex_t IT_0819 = IT_0301*IT_0818;
    const ccomplex_t IT_0820 = IT_0329 + IT_0535 + IT_0630 + IT_0634;
    const ccomplex_t IT_0821 = (-2)*IT_0338 + -IT_0534 + -IT_0628 + -IT_0633;
    const ccomplex_t IT_0822 = IT_0820 + IT_0821;
    const ccomplex_t IT_0823 = IT_0326*IT_0822;
    const ccomplex_t IT_0824 = IT_0350 + IT_0354 + IT_0357 + IT_0360 + IT_0361
       + IT_0362;
    const ccomplex_t IT_0825 = IT_0349*IT_0824;
    const ccomplex_t IT_0826 = IT_0371 + IT_0375 + IT_0377 + IT_0381 + IT_0382;
    const ccomplex_t IT_0827 = IT_0370*IT_0826;
    const ccomplex_t IT_0828 = IT_0350 + IT_0354 + IT_0357 + IT_0360 + IT_0361;
    const ccomplex_t IT_0829 = IT_0390*IT_0828;
    const ccomplex_t IT_0830 = IT_0395*IT_0826;
    const ccomplex_t IT_0831 = IT_0398*IT_0828;
    const ccomplex_t IT_0832 = IT_0403 + IT_0405 + IT_0411 + IT_0412 + IT_0413
       + IT_0414;
    const ccomplex_t IT_0833 = IT_0401*IT_0832;
    const ccomplex_t IT_0834 = IT_0420 + IT_0424 + IT_0426 + IT_0428 + IT_0429;
    const ccomplex_t IT_0835 = IT_0419*IT_0834;
    const ccomplex_t IT_0836 = IT_0403 + IT_0405 + IT_0412 + IT_0413 + IT_0414;
    const ccomplex_t IT_0837 = IT_0436*IT_0836;
    const ccomplex_t IT_0838 = IT_0441*IT_0834;
    const ccomplex_t IT_0839 = IT_0444*IT_0836;
    const ccomplex_t IT_0840 = IT_0449 + IT_0451 + IT_0457 + IT_0458 + IT_0459
       + IT_0462;
    const ccomplex_t IT_0841 = IT_0447*IT_0840;
    const ccomplex_t IT_0842 = IT_0475 + IT_0477 + IT_0481 + IT_0482 + IT_0483;
    const ccomplex_t IT_0843 = IT_0474*IT_0842;
    const ccomplex_t IT_0844 = IT_0449 + IT_0451 + IT_0457 + IT_0458 + IT_0459;
    const ccomplex_t IT_0845 = IT_0496*IT_0844;
    const ccomplex_t IT_0846 = IT_0501*IT_0842;
    const ccomplex_t IT_0847 = IT_0507*IT_0844;
    const ccomplex_t IT_0848 = IT_0513*IT_0798;
    const ccomplex_t IT_0849 = IT_0522*IT_0802;
    const ccomplex_t IT_0850 = IT_0531*IT_0806;
    const ccomplex_t IT_0851 = -IT_0563 + -IT_0564 + -IT_0565 + -IT_0568 + 
      -IT_0569 + -IT_0570 + 2*IT_0572 + (-2)*IT_0573 + (-2)*IT_0574 + (-2)
      *IT_0577 + (-2)*IT_0580 + (-2)*IT_0583 + (-2)*IT_0584 + (-2)*IT_0585 + 2
      *IT_0586 + 2*IT_0642 + 2*IT_0647 + (-2)*IT_0652 + 2*IT_0662 + (-2)*IT_0668
       + 2*IT_0673 + -IT_0809 + -IT_0810 + -IT_0811 + 2*IT_0815 + 2*IT_0819 + 2
      *IT_0823 + 4*IT_0825 + 2*IT_0827 + 2*IT_0829 + (-2)*IT_0830 + (-2)*IT_0831
       + 4*IT_0833 + 2*IT_0835 + 2*IT_0837 + (-2)*IT_0838 + (-2)*IT_0839 + 4
      *IT_0841 + 2*IT_0843 + 2*IT_0845 + (-2)*IT_0846 + (-2)*IT_0847 + -IT_0848 
      + -IT_0849 + -IT_0850;
    const ccomplex_t IT_0852 = IT_0808 + IT_0851;
    const ccomplex_t IT_0853 = IT_0005*IT_0852;
    const ccomplex_t IT_0854 = 1.4142135623731*IT_0853;
    const ccomplex_t IT_0855 = 0.125*IT_0854;
    return create_ccomplex_return((0 + _Complex_I*0.25)*IT_0552 + (0 +
       _Complex_I*0.25)*IT_0696 + (0 + _Complex_I*0.25)*IT_0783 + (0 +
       _Complex_I*0.25)*IT_0855);
}

