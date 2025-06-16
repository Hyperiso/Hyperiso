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
    const creal_t e_em = param->e_em;
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
    const ccomplex_t IT_0006 = pow(m_b, 2);
    const ccomplex_t IT_0007 = pow(m_s, 2);
    const ccomplex_t IT_0008 = pow(m_u, 2);
    const ccomplex_t IT_0009 = pow(M_W, 2);
    const ccomplex_t IT_0010 = mtylt_C0iC(3, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0011 = mtylt_C0iC(6, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0012 = mtylt_C0iC(0, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0013 = IT_0010 + IT_0011 + IT_0012;
    const ccomplex_t IT_0014 = 0.101321183642338*m_u;
    const ccomplex_t IT_0015 = pow(M_W, -1);
    const ccomplex_t IT_0016 = sin(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0021 = 0.666666666666667*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_u*V_ub*e_em
      *IT_0015*IT_0017;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = IT_0019*IT_0021*IT_0023;
    const ccomplex_t IT_0025 = IT_0014*IT_0024;
    const ccomplex_t IT_0026 = IT_0013*IT_0025;
    const ccomplex_t IT_0027 = pow(m_c, 2);
    const ccomplex_t IT_0028 = mtylt_C0iC(3, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0027, IT_0009, IT_0027, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0029 = mtylt_C0iC(6, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0027, IT_0009, IT_0027, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0030 = mtylt_C0iC(0, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0027, IT_0009, IT_0027, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0031 = IT_0028 + IT_0029 + IT_0030;
    const ccomplex_t IT_0032 = 0.101321183642338*m_c;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0015*IT_0017;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0015*IT_0017;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0021*IT_0034*IT_0036;
    const ccomplex_t IT_0038 = IT_0032*IT_0037;
    const ccomplex_t IT_0039 = IT_0031*IT_0038;
    const ccomplex_t IT_0040 = pow(m_t, 2);
    const ccomplex_t IT_0041 = mtylt_C0iC(3, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0040, IT_0009, IT_0040, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0042 = mtylt_C0iC(6, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0040, IT_0009, IT_0040, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0043 = mtylt_C0iC(0, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0040, IT_0009, IT_0040, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0044 = IT_0041 + IT_0042 + IT_0043;
    const ccomplex_t IT_0045 = 0.101321183642338*m_t;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0015*IT_0017;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0015*IT_0017;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = IT_0021*IT_0047*IT_0049;
    const ccomplex_t IT_0051 = IT_0045*IT_0050;
    const ccomplex_t IT_0052 = IT_0044*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0015*IT_0017;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0015*IT_0017;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0021*IT_0054*IT_0056;
    const ccomplex_t IT_0058 = IT_0045*IT_0057;
    const ccomplex_t IT_0059 = IT_0044*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0017;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*V_ub*e_em
      *IT_0017;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = IT_0061*IT_0063;
    const ccomplex_t IT_0065 = 0.101321183642338*IT_0064;
    const ccomplex_t IT_0066 = mtylt_C0iC(15, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0067 = (-2)*IT_0020;
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = m_s*IT_0068;
    const ccomplex_t IT_0070 = mtylt_C0iC(6, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0071 = -IT_0020;
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = m_b*IT_0072;
    const ccomplex_t IT_0074 = mtylt_C0iC(18, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0075 = IT_0067*IT_0074;
    const ccomplex_t IT_0076 = m_b*IT_0075;
    const ccomplex_t IT_0077 = mtylt_C0iC(3, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0078 = IT_0071*IT_0077;
    const ccomplex_t IT_0079 = m_s*IT_0078;
    const ccomplex_t IT_0080 = 2*IT_0020;
    const ccomplex_t IT_0081 = IT_0077*IT_0080;
    const ccomplex_t IT_0082 = m_s*IT_0081;
    const ccomplex_t IT_0083 = IT_0070*IT_0080;
    const ccomplex_t IT_0084 = m_s*IT_0083;
    const ccomplex_t IT_0085 = IT_0069 + IT_0073 + IT_0076 + IT_0079 + IT_0082
       + IT_0084;
    const ccomplex_t IT_0086 = IT_0020*IT_0070;
    const ccomplex_t IT_0087 = m_s*IT_0086;
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0085 + IT_0088;
    const ccomplex_t IT_0090 = IT_0065*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0017;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0017;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0092*IT_0094;
    const ccomplex_t IT_0096 = 0.101321183642338*IT_0095;
    const ccomplex_t IT_0097 = mtylt_C0iC(3, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0027, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0098 = IT_0071*IT_0097;
    const ccomplex_t IT_0099 = m_s*IT_0098;
    const ccomplex_t IT_0100 = IT_0080*IT_0097;
    const ccomplex_t IT_0101 = m_s*IT_0100;
    const ccomplex_t IT_0102 = mtylt_C0iC(18, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0027, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0103 = IT_0067*IT_0102;
    const ccomplex_t IT_0104 = m_b*IT_0103;
    const ccomplex_t IT_0105 = mtylt_C0iC(15, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0027, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0106 = IT_0067*IT_0105;
    const ccomplex_t IT_0107 = m_s*IT_0106;
    const ccomplex_t IT_0108 = mtylt_C0iC(6, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0027, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0109 = IT_0080*IT_0108;
    const ccomplex_t IT_0110 = m_s*IT_0109;
    const ccomplex_t IT_0111 = IT_0071*IT_0108;
    const ccomplex_t IT_0112 = m_b*IT_0111;
    const ccomplex_t IT_0113 = IT_0099 + IT_0101 + IT_0104 + IT_0107 + IT_0110
       + IT_0112;
    const ccomplex_t IT_0114 = IT_0020*IT_0108;
    const ccomplex_t IT_0115 = m_s*IT_0114;
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = IT_0113 + IT_0116;
    const ccomplex_t IT_0118 = IT_0096*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0017;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0017;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0120*IT_0122;
    const ccomplex_t IT_0124 = 0.101321183642338*IT_0123;
    const ccomplex_t IT_0125 = mtylt_C0iC(15, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0040, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0126 = IT_0067*IT_0125;
    const ccomplex_t IT_0127 = m_s*IT_0126;
    const ccomplex_t IT_0128 = mtylt_C0iC(18, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0040, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0129 = IT_0067*IT_0128;
    const ccomplex_t IT_0130 = m_b*IT_0129;
    const ccomplex_t IT_0131 = mtylt_C0iC(3, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0040, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0132 = IT_0071*IT_0131;
    const ccomplex_t IT_0133 = m_s*IT_0132;
    const ccomplex_t IT_0134 = mtylt_C0iC(6, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0040, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0135 = IT_0071*IT_0134;
    const ccomplex_t IT_0136 = m_b*IT_0135;
    const ccomplex_t IT_0137 = IT_0080*IT_0131;
    const ccomplex_t IT_0138 = m_s*IT_0137;
    const ccomplex_t IT_0139 = IT_0080*IT_0134;
    const ccomplex_t IT_0140 = m_s*IT_0139;
    const ccomplex_t IT_0141 = IT_0127 + IT_0130 + IT_0133 + IT_0136 + IT_0138
       + IT_0140;
    const ccomplex_t IT_0142 = IT_0020*IT_0134;
    const ccomplex_t IT_0143 = m_s*IT_0142;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0141 + IT_0144;
    const ccomplex_t IT_0146 = IT_0124*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0015*IT_0017;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = IT_0021*IT_0023*IT_0148;
    const ccomplex_t IT_0150 = 0.101321183642338*IT_0149;
    const ccomplex_t IT_0151 = m_s*IT_0011;
    const ccomplex_t IT_0152 = mtylt_C0iC(15, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0153 = m_s*IT_0152;
    const ccomplex_t IT_0154 = mtylt_C0iC(18, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0155 = m_s*IT_0154;
    const ccomplex_t IT_0156 = IT_0151 + IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = m_b*IT_0011;
    const ccomplex_t IT_0158 = mtylt_C0iC(12, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0159 = m_b*IT_0158;
    const ccomplex_t IT_0160 = m_b*IT_0152;
    const ccomplex_t IT_0161 = m_b*IT_0154;
    const ccomplex_t IT_0162 = m_b*IT_0010;
    const ccomplex_t IT_0163 = -IT_0157 + -IT_0159 + (-2)*IT_0160 + -IT_0161 +
       -IT_0162;
    const ccomplex_t IT_0164 = IT_0156 + IT_0163;
    const ccomplex_t IT_0165 = IT_0150*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1.4142135623731)*m_b*V_ub*e_em
      *IT_0015*IT_0017;
    const ccomplex_t IT_0167 = (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = IT_0019*IT_0021*IT_0167;
    const ccomplex_t IT_0169 = 0.101321183642338*IT_0168;
    const ccomplex_t IT_0170 = IT_0164*IT_0169;
    const ccomplex_t IT_0171 = IT_0021*IT_0148*IT_0167;
    const ccomplex_t IT_0172 = IT_0014*IT_0171;
    const ccomplex_t IT_0173 = IT_0013*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0015*IT_0017;
    const ccomplex_t IT_0175 = 0.5*IT_0174;
    const ccomplex_t IT_0176 = IT_0021*IT_0036*IT_0175;
    const ccomplex_t IT_0177 = 0.101321183642338*IT_0176;
    const ccomplex_t IT_0178 = m_s*IT_0029;
    const ccomplex_t IT_0179 = mtylt_C0iC(15, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0027, IT_0009, IT_0027, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0180 = m_s*IT_0179;
    const ccomplex_t IT_0181 = mtylt_C0iC(18, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0027, IT_0009, IT_0027, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0182 = m_s*IT_0181;
    const ccomplex_t IT_0183 = IT_0178 + IT_0180 + IT_0182;
    const ccomplex_t IT_0184 = m_b*IT_0029;
    const ccomplex_t IT_0185 = mtylt_C0iC(12, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0027, IT_0009, IT_0027, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0186 = m_b*IT_0185;
    const ccomplex_t IT_0187 = m_b*IT_0179;
    const ccomplex_t IT_0188 = m_b*IT_0181;
    const ccomplex_t IT_0189 = m_b*IT_0028;
    const ccomplex_t IT_0190 = -IT_0184 + -IT_0186 + (-2)*IT_0187 + -IT_0188 +
       -IT_0189;
    const ccomplex_t IT_0191 = IT_0183 + IT_0190;
    const ccomplex_t IT_0192 = IT_0177*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0015*IT_0017;
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = IT_0021*IT_0034*IT_0194;
    const ccomplex_t IT_0196 = 0.101321183642338*IT_0195;
    const ccomplex_t IT_0197 = IT_0191*IT_0196;
    const ccomplex_t IT_0198 = IT_0021*IT_0175*IT_0194;
    const ccomplex_t IT_0199 = IT_0032*IT_0198;
    const ccomplex_t IT_0200 = IT_0031*IT_0199;
    const ccomplex_t IT_0201 = IT_0021*IT_0049*IT_0056;
    const ccomplex_t IT_0202 = 0.101321183642338*IT_0201;
    const ccomplex_t IT_0203 = m_s*IT_0042;
    const ccomplex_t IT_0204 = mtylt_C0iC(15, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0040, IT_0009, IT_0040, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0205 = m_s*IT_0204;
    const ccomplex_t IT_0206 = mtylt_C0iC(18, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0040, IT_0009, IT_0040, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0207 = m_s*IT_0206;
    const ccomplex_t IT_0208 = IT_0203 + IT_0205 + IT_0207;
    const ccomplex_t IT_0209 = mtylt_C0iC(12, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0040, IT_0009, IT_0040, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0210 = m_b*IT_0209;
    const ccomplex_t IT_0211 = m_b*IT_0204;
    const ccomplex_t IT_0212 = m_b*IT_0206;
    const ccomplex_t IT_0213 = m_b*IT_0041;
    const ccomplex_t IT_0214 = m_b*IT_0042;
    const ccomplex_t IT_0215 = -IT_0210 + (-2)*IT_0211 + -IT_0212 + -IT_0213 +
       -IT_0214;
    const ccomplex_t IT_0216 = IT_0208 + IT_0215;
    const ccomplex_t IT_0217 = IT_0202*IT_0216;
    const ccomplex_t IT_0218 = IT_0021*IT_0047*IT_0054;
    const ccomplex_t IT_0219 = 0.101321183642338*IT_0218;
    const ccomplex_t IT_0220 = IT_0216*IT_0219;
    const ccomplex_t IT_0221 = IT_0026 + IT_0039 + IT_0052 + IT_0059 + IT_0090
       + IT_0118 + IT_0146 + IT_0165 + IT_0170 + IT_0173 + IT_0192 + IT_0197 +
       IT_0200 + IT_0217 + IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*M_W*e_em;
    const ccomplex_t IT_0223 = IT_0034*IT_0094*IT_0222;
    const ccomplex_t IT_0224 = 0.101321183642338*IT_0223;
    const ccomplex_t IT_0225 = IT_0108*IT_0224;
    const ccomplex_t IT_0226 = IT_0047*IT_0122*IT_0222;
    const ccomplex_t IT_0227 = 0.101321183642338*IT_0226;
    const ccomplex_t IT_0228 = IT_0134*IT_0227;
    const ccomplex_t IT_0229 = IT_0034*IT_0036;
    const ccomplex_t IT_0230 = IT_0032*IT_0229;
    const ccomplex_t IT_0231 = mtylt_C0iC(0, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0027, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0232 = IT_0071*IT_0231;
    const ccomplex_t IT_0233 = IT_0067*IT_0231;
    const ccomplex_t IT_0234 = IT_0067*IT_0108;
    const ccomplex_t IT_0235 = -IT_0233 + -IT_0234;
    const ccomplex_t IT_0236 = IT_0232 + IT_0235;
    const ccomplex_t IT_0237 = IT_0230*IT_0236;
    const ccomplex_t IT_0238 = IT_0148*IT_0167;
    const ccomplex_t IT_0239 = IT_0014*IT_0238;
    const ccomplex_t IT_0240 = mtylt_C0iC(0, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0241 = IT_0071*IT_0240;
    const ccomplex_t IT_0242 = IT_0067*IT_0240;
    const ccomplex_t IT_0243 = IT_0067*IT_0070;
    const ccomplex_t IT_0244 = -IT_0242 + -IT_0243;
    const ccomplex_t IT_0245 = IT_0241 + IT_0244;
    const ccomplex_t IT_0246 = IT_0239*IT_0245;
    const ccomplex_t IT_0247 = IT_0019*IT_0063*IT_0222;
    const ccomplex_t IT_0248 = 0.101321183642338*IT_0247;
    const ccomplex_t IT_0249 = IT_0070*IT_0248;
    const ccomplex_t IT_0250 = IT_0021*IT_0061*IT_0063;
    const ccomplex_t IT_0251 = 0.101321183642338*IT_0250;
    const ccomplex_t IT_0252 = m_s*IT_0010;
    const ccomplex_t IT_0253 = IT_0151 + IT_0153 + IT_0155 + IT_0252;
    const ccomplex_t IT_0254 = IT_0163 + IT_0253;
    const ccomplex_t IT_0255 = IT_0251*IT_0254;
    const ccomplex_t IT_0256 = IT_0021*IT_0092*IT_0094;
    const ccomplex_t IT_0257 = 0.101321183642338*IT_0256;
    const ccomplex_t IT_0258 = m_s*IT_0028;
    const ccomplex_t IT_0259 = IT_0178 + IT_0180 + IT_0182 + IT_0258;
    const ccomplex_t IT_0260 = IT_0190 + IT_0259;
    const ccomplex_t IT_0261 = IT_0257*IT_0260;
    const ccomplex_t IT_0262 = IT_0021*IT_0120*IT_0122;
    const ccomplex_t IT_0263 = 0.101321183642338*IT_0262;
    const ccomplex_t IT_0264 = m_s*IT_0041;
    const ccomplex_t IT_0265 = IT_0203 + IT_0205 + IT_0207 + IT_0264;
    const ccomplex_t IT_0266 = IT_0215 + IT_0265;
    const ccomplex_t IT_0267 = IT_0263*IT_0266;
    const ccomplex_t IT_0268 = IT_0019*IT_0023;
    const ccomplex_t IT_0269 = IT_0014*IT_0268;
    const ccomplex_t IT_0270 = IT_0245*IT_0269;
    const ccomplex_t IT_0271 = IT_0023*IT_0148;
    const ccomplex_t IT_0272 = 0.101321183642338*IT_0271;
    const ccomplex_t IT_0273 = IT_0073 + IT_0079;
    const ccomplex_t IT_0274 = IT_0067*IT_0077;
    const ccomplex_t IT_0275 = m_s*IT_0274;
    const ccomplex_t IT_0276 = m_b*IT_0243;
    const ccomplex_t IT_0277 = -IT_0069 + -IT_0076 + -IT_0275 + -IT_0276;
    const ccomplex_t IT_0278 = IT_0273 + IT_0277;
    const ccomplex_t IT_0279 = IT_0272*IT_0278;
    const ccomplex_t IT_0280 = IT_0036*IT_0175;
    const ccomplex_t IT_0281 = 0.101321183642338*IT_0280;
    const ccomplex_t IT_0282 = IT_0099 + IT_0112;
    const ccomplex_t IT_0283 = m_b*IT_0234;
    const ccomplex_t IT_0284 = IT_0067*IT_0097;
    const ccomplex_t IT_0285 = m_s*IT_0284;
    const ccomplex_t IT_0286 = -IT_0104 + -IT_0107 + -IT_0283 + -IT_0285;
    const ccomplex_t IT_0287 = IT_0282 + IT_0286;
    const ccomplex_t IT_0288 = IT_0281*IT_0287;
    const ccomplex_t IT_0289 = IT_0047*IT_0049;
    const ccomplex_t IT_0290 = IT_0045*IT_0289;
    const ccomplex_t IT_0291 = mtylt_C0iC(0, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0040, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0292 = IT_0071*IT_0291;
    const ccomplex_t IT_0293 = IT_0067*IT_0291;
    const ccomplex_t IT_0294 = IT_0067*IT_0134;
    const ccomplex_t IT_0295 = -IT_0293 + -IT_0294;
    const ccomplex_t IT_0296 = IT_0292 + IT_0295;
    const ccomplex_t IT_0297 = IT_0290*IT_0296;
    const ccomplex_t IT_0298 = IT_0049*IT_0056;
    const ccomplex_t IT_0299 = 0.101321183642338*IT_0298;
    const ccomplex_t IT_0300 = IT_0133 + IT_0136;
    const ccomplex_t IT_0301 = IT_0067*IT_0131;
    const ccomplex_t IT_0302 = m_s*IT_0301;
    const ccomplex_t IT_0303 = m_b*IT_0294;
    const ccomplex_t IT_0304 = -IT_0127 + -IT_0130 + -IT_0302 + -IT_0303;
    const ccomplex_t IT_0305 = IT_0300 + IT_0304;
    const ccomplex_t IT_0306 = IT_0299*IT_0305;
    const ccomplex_t IT_0307 = IT_0019*IT_0167;
    const ccomplex_t IT_0308 = 0.101321183642338*IT_0307;
    const ccomplex_t IT_0309 = IT_0278*IT_0308;
    const ccomplex_t IT_0310 = IT_0034*IT_0194;
    const ccomplex_t IT_0311 = 0.101321183642338*IT_0310;
    const ccomplex_t IT_0312 = IT_0287*IT_0311;
    const ccomplex_t IT_0313 = IT_0175*IT_0194;
    const ccomplex_t IT_0314 = IT_0032*IT_0313;
    const ccomplex_t IT_0315 = IT_0236*IT_0314;
    const ccomplex_t IT_0316 = IT_0047*IT_0054;
    const ccomplex_t IT_0317 = 0.101321183642338*IT_0316;
    const ccomplex_t IT_0318 = IT_0305*IT_0317;
    const ccomplex_t IT_0319 = IT_0054*IT_0056;
    const ccomplex_t IT_0320 = IT_0045*IT_0319;
    const ccomplex_t IT_0321 = IT_0296*IT_0320;
    const ccomplex_t IT_0322 = -IT_0225 + -IT_0228 + 0.5*IT_0237 + 0.5*IT_0246
       + -IT_0249 + 2*IT_0255 + 2*IT_0261 + 2*IT_0267 + 0.5*IT_0270 + (-0.5)
      *IT_0279 + (-0.5)*IT_0288 + 0.5*IT_0297 + (-0.5)*IT_0306 + (-0.5)*IT_0309 
      + (-0.5)*IT_0312 + 0.5*IT_0315 + (-0.5)*IT_0318 + 0.5*IT_0321;
    const ccomplex_t IT_0323 = IT_0221 + IT_0322;
    const ccomplex_t IT_0324 = IT_0005*IT_0323;
    const ccomplex_t IT_0325 = 1.4142135623731*IT_0324;
    const ccomplex_t IT_0326 = 0.25*IT_0325;
    const ccomplex_t IT_0327 = IT_0203 + IT_0207;
    const ccomplex_t IT_0328 = -IT_0211 + -IT_0212 + -IT_0214;
    const ccomplex_t IT_0329 = IT_0327 + IT_0328;
    const ccomplex_t IT_0330 = IT_0202*IT_0329;
    const ccomplex_t IT_0331 = IT_0219*IT_0329;
    const ccomplex_t IT_0332 = IT_0061*IT_0167*IT_0222;
    const ccomplex_t IT_0333 = 0.101321183642338*IT_0332;
    const ccomplex_t IT_0334 = IT_0077*IT_0333;
    const ccomplex_t IT_0335 = IT_0092*IT_0194*IT_0222;
    const ccomplex_t IT_0336 = 0.101321183642338*IT_0335;
    const ccomplex_t IT_0337 = IT_0097*IT_0336;
    const ccomplex_t IT_0338 = IT_0054*IT_0120*IT_0222;
    const ccomplex_t IT_0339 = 0.101321183642338*IT_0338;
    const ccomplex_t IT_0340 = IT_0131*IT_0339;
    const ccomplex_t IT_0341 = IT_0042*IT_0058;
    const ccomplex_t IT_0342 = IT_0011*IT_0025;
    const ccomplex_t IT_0343 = IT_0011*IT_0172;
    const ccomplex_t IT_0344 = IT_0029*IT_0038;
    const ccomplex_t IT_0345 = IT_0029*IT_0199;
    const ccomplex_t IT_0346 = IT_0042*IT_0051;
    const ccomplex_t IT_0347 = m_b*IT_0078;
    const ccomplex_t IT_0348 = IT_0079 + IT_0347;
    const ccomplex_t IT_0349 = mtylt_C0iC(12, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0350 = IT_0067*IT_0349;
    const ccomplex_t IT_0351 = m_s*IT_0350;
    const ccomplex_t IT_0352 = m_b*IT_0068;
    const ccomplex_t IT_0353 = m_b*IT_0086;
    const ccomplex_t IT_0354 = 2*IT_0073 + -IT_0275 + -IT_0276 + -IT_0351 + 
      -IT_0352 + -IT_0353;
    const ccomplex_t IT_0355 = IT_0348 + IT_0354;
    const ccomplex_t IT_0356 = IT_0065*IT_0355;
    const ccomplex_t IT_0357 = m_b*IT_0098;
    const ccomplex_t IT_0358 = IT_0099 + IT_0357;
    const ccomplex_t IT_0359 = mtylt_C0iC(12, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0027, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0360 = IT_0067*IT_0359;
    const ccomplex_t IT_0361 = m_s*IT_0360;
    const ccomplex_t IT_0362 = m_b*IT_0106;
    const ccomplex_t IT_0363 = m_b*IT_0114;
    const ccomplex_t IT_0364 = 2*IT_0112 + -IT_0283 + -IT_0285 + -IT_0361 + 
      -IT_0362 + -IT_0363;
    const ccomplex_t IT_0365 = IT_0358 + IT_0364;
    const ccomplex_t IT_0366 = IT_0096*IT_0365;
    const ccomplex_t IT_0367 = m_b*IT_0132;
    const ccomplex_t IT_0368 = IT_0133 + IT_0367;
    const ccomplex_t IT_0369 = mtylt_C0iC(12, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0040, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0370 = IT_0067*IT_0369;
    const ccomplex_t IT_0371 = m_s*IT_0370;
    const ccomplex_t IT_0372 = m_b*IT_0142;
    const ccomplex_t IT_0373 = m_b*IT_0126;
    const ccomplex_t IT_0374 = 2*IT_0136 + -IT_0302 + -IT_0303 + -IT_0371 + 
      -IT_0372 + -IT_0373;
    const ccomplex_t IT_0375 = IT_0368 + IT_0374;
    const ccomplex_t IT_0376 = IT_0124*IT_0375;
    const ccomplex_t IT_0377 = IT_0151 + IT_0155;
    const ccomplex_t IT_0378 = -IT_0157 + -IT_0160 + -IT_0161;
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = IT_0150*IT_0379;
    const ccomplex_t IT_0381 = IT_0169*IT_0379;
    const ccomplex_t IT_0382 = IT_0178 + IT_0182;
    const ccomplex_t IT_0383 = -IT_0184 + -IT_0187 + -IT_0188;
    const ccomplex_t IT_0384 = IT_0382 + IT_0383;
    const ccomplex_t IT_0385 = IT_0177*IT_0384;
    const ccomplex_t IT_0386 = IT_0196*IT_0384;
    const ccomplex_t IT_0387 = IT_0330 + IT_0331 + IT_0334 + IT_0337 + IT_0340
       + IT_0341 + IT_0342 + IT_0343 + IT_0344 + IT_0345 + IT_0346 + IT_0356 +
       IT_0366 + IT_0376 + IT_0380 + IT_0381 + IT_0385 + IT_0386;
    const ccomplex_t IT_0388 = -IT_0157 + -IT_0160 + -IT_0161 + -IT_0162;
    const ccomplex_t IT_0389 = IT_0377 + IT_0388;
    const ccomplex_t IT_0390 = IT_0251*IT_0389;
    const ccomplex_t IT_0391 = -IT_0184 + -IT_0187 + -IT_0188 + -IT_0189;
    const ccomplex_t IT_0392 = IT_0382 + IT_0391;
    const ccomplex_t IT_0393 = IT_0257*IT_0392;
    const ccomplex_t IT_0394 = -IT_0211 + -IT_0212 + -IT_0213 + -IT_0214;
    const ccomplex_t IT_0395 = IT_0327 + IT_0394;
    const ccomplex_t IT_0396 = IT_0263*IT_0395;
    const ccomplex_t IT_0397 = IT_0241 + IT_0274;
    const ccomplex_t IT_0398 = IT_0269*IT_0397;
    const ccomplex_t IT_0399 = IT_0073 + IT_0079 + IT_0351 + IT_0352;
    const ccomplex_t IT_0400 = IT_0272*IT_0399;
    const ccomplex_t IT_0401 = IT_0232 + IT_0284;
    const ccomplex_t IT_0402 = IT_0230*IT_0401;
    const ccomplex_t IT_0403 = IT_0099 + IT_0112 + IT_0361 + IT_0362;
    const ccomplex_t IT_0404 = IT_0281*IT_0403;
    const ccomplex_t IT_0405 = IT_0292 + IT_0301;
    const ccomplex_t IT_0406 = IT_0290*IT_0405;
    const ccomplex_t IT_0407 = IT_0133 + IT_0136 + IT_0371 + IT_0373;
    const ccomplex_t IT_0408 = IT_0299*IT_0407;
    const ccomplex_t IT_0409 = IT_0308*IT_0399;
    const ccomplex_t IT_0410 = IT_0239*IT_0397;
    const ccomplex_t IT_0411 = IT_0311*IT_0403;
    const ccomplex_t IT_0412 = IT_0314*IT_0401;
    const ccomplex_t IT_0413 = IT_0317*IT_0407;
    const ccomplex_t IT_0414 = IT_0320*IT_0405;
    const ccomplex_t IT_0415 = 2*IT_0390 + 2*IT_0393 + 2*IT_0396 + 0.5*IT_0398
       + (-0.5)*IT_0400 + 0.5*IT_0402 + (-0.5)*IT_0404 + 0.5*IT_0406 + (-0.5)
      *IT_0408 + (-0.5)*IT_0409 + 0.5*IT_0410 + (-0.5)*IT_0411 + 0.5*IT_0412 + (
      -0.5)*IT_0413 + 0.5*IT_0414;
    const ccomplex_t IT_0416 = IT_0387 + IT_0415;
    const ccomplex_t IT_0417 = IT_0005*IT_0416;
    const ccomplex_t IT_0418 = 1.4142135623731*IT_0417;
    const ccomplex_t IT_0419 = (-0.25)*IT_0418;
    const ccomplex_t IT_0420 = 2*IT_0160;
    const ccomplex_t IT_0421 = IT_0151 + IT_0153 + IT_0155 + IT_0157 + IT_0159
       + IT_0161 + IT_0162;
    const ccomplex_t IT_0422 = IT_0420 + IT_0421;
    const ccomplex_t IT_0423 = IT_0150*IT_0422;
    const ccomplex_t IT_0424 = 2*IT_0187;
    const ccomplex_t IT_0425 = IT_0178 + IT_0180 + IT_0182 + IT_0184 + IT_0186
       + IT_0188 + IT_0189;
    const ccomplex_t IT_0426 = IT_0424 + IT_0425;
    const ccomplex_t IT_0427 = IT_0177*IT_0426;
    const ccomplex_t IT_0428 = 2*IT_0211;
    const ccomplex_t IT_0429 = IT_0203 + IT_0205 + IT_0207 + IT_0210 + IT_0212
       + IT_0213 + IT_0214;
    const ccomplex_t IT_0430 = IT_0428 + IT_0429;
    const ccomplex_t IT_0431 = IT_0202*IT_0430;
    const ccomplex_t IT_0432 = IT_0069 + IT_0079 + IT_0082 + IT_0084;
    const ccomplex_t IT_0433 = -IT_0073 + -IT_0076 + -IT_0087;
    const ccomplex_t IT_0434 = IT_0432 + IT_0433;
    const ccomplex_t IT_0435 = IT_0065*IT_0434;
    const ccomplex_t IT_0436 = IT_0099 + IT_0101 + IT_0107 + IT_0110;
    const ccomplex_t IT_0437 = -IT_0104 + -IT_0112 + -IT_0115;
    const ccomplex_t IT_0438 = IT_0436 + IT_0437;
    const ccomplex_t IT_0439 = IT_0096*IT_0438;
    const ccomplex_t IT_0440 = IT_0127 + IT_0133 + IT_0138 + IT_0140;
    const ccomplex_t IT_0441 = -IT_0130 + -IT_0136 + -IT_0143;
    const ccomplex_t IT_0442 = IT_0440 + IT_0441;
    const ccomplex_t IT_0443 = IT_0124*IT_0442;
    const ccomplex_t IT_0444 = IT_0026 + IT_0039 + IT_0052 + IT_0423 + IT_0427
       + IT_0431 + IT_0435 + IT_0439 + IT_0443;
    const ccomplex_t IT_0445 = IT_0151 + IT_0153 + IT_0155 + IT_0157 + IT_0159
       + IT_0161 + IT_0162 + IT_0252;
    const ccomplex_t IT_0446 = IT_0420 + IT_0445;
    const ccomplex_t IT_0447 = IT_0251*IT_0446;
    const ccomplex_t IT_0448 = IT_0169*IT_0422;
    const ccomplex_t IT_0449 = IT_0178 + IT_0180 + IT_0182 + IT_0184 + IT_0186
       + IT_0188 + IT_0189 + IT_0258;
    const ccomplex_t IT_0450 = IT_0424 + IT_0449;
    const ccomplex_t IT_0451 = IT_0257*IT_0450;
    const ccomplex_t IT_0452 = IT_0196*IT_0426;
    const ccomplex_t IT_0453 = IT_0203 + IT_0205 + IT_0207 + IT_0210 + IT_0212
       + IT_0213 + IT_0214 + IT_0264;
    const ccomplex_t IT_0454 = IT_0428 + IT_0453;
    const ccomplex_t IT_0455 = IT_0263*IT_0454;
    const ccomplex_t IT_0456 = IT_0219*IT_0430;
    const ccomplex_t IT_0457 = IT_0076 + IT_0079 + IT_0276;
    const ccomplex_t IT_0458 = -IT_0069 + -IT_0073 + -IT_0275;
    const ccomplex_t IT_0459 = IT_0457 + IT_0458;
    const ccomplex_t IT_0460 = IT_0272*IT_0459;
    const ccomplex_t IT_0461 = -IT_0107 + -IT_0112 + -IT_0285;
    const ccomplex_t IT_0462 = IT_0099 + IT_0104 + IT_0283;
    const ccomplex_t IT_0463 = IT_0461 + IT_0462;
    const ccomplex_t IT_0464 = IT_0281*IT_0463;
    const ccomplex_t IT_0465 = IT_0130 + IT_0133 + IT_0303;
    const ccomplex_t IT_0466 = -IT_0127 + -IT_0136 + -IT_0302;
    const ccomplex_t IT_0467 = IT_0465 + IT_0466;
    const ccomplex_t IT_0468 = IT_0299*IT_0467;
    const ccomplex_t IT_0469 = IT_0308*IT_0459;
    const ccomplex_t IT_0470 = IT_0311*IT_0463;
    const ccomplex_t IT_0471 = IT_0317*IT_0467;
    const ccomplex_t IT_0472 = -IT_0059 + -IT_0173 + -IT_0200 + -IT_0225 + 
      -IT_0228 + 0.5*IT_0237 + (-0.5)*IT_0246 + -IT_0249 + 0.5*IT_0270 + 0.5
      *IT_0297 + (-0.5)*IT_0315 + (-0.5)*IT_0321 + 2*IT_0447 + -IT_0448 + 2
      *IT_0451 + -IT_0452 + 2*IT_0455 + -IT_0456 + (-0.5)*IT_0460 + (-0.5)
      *IT_0464 + (-0.5)*IT_0468 + 0.5*IT_0469 + 0.5*IT_0470 + 0.5*IT_0471;
    const ccomplex_t IT_0473 = IT_0444 + IT_0472;
    const ccomplex_t IT_0474 = IT_0005*IT_0473;
    const ccomplex_t IT_0475 = 1.4142135623731*IT_0474;
    const ccomplex_t IT_0476 = (-0.25)*IT_0475;
    const ccomplex_t IT_0477 = IT_0151 + IT_0155 + IT_0157 + IT_0160 + IT_0161;
    const ccomplex_t IT_0478 = IT_0150*IT_0477;
    const ccomplex_t IT_0479 = IT_0203 + IT_0207 + IT_0211 + IT_0212 + IT_0214;
    const ccomplex_t IT_0480 = IT_0202*IT_0479;
    const ccomplex_t IT_0481 = IT_0079 + IT_0276 + IT_0352 + IT_0353;
    const ccomplex_t IT_0482 = (-2)*IT_0073 + -IT_0275 + -IT_0347 + -IT_0351;
    const ccomplex_t IT_0483 = IT_0481 + IT_0482;
    const ccomplex_t IT_0484 = IT_0065*IT_0483;
    const ccomplex_t IT_0485 = IT_0099 + IT_0283 + IT_0362 + IT_0363;
    const ccomplex_t IT_0486 = (-2)*IT_0112 + -IT_0285 + -IT_0357 + -IT_0361;
    const ccomplex_t IT_0487 = IT_0485 + IT_0486;
    const ccomplex_t IT_0488 = IT_0096*IT_0487;
    const ccomplex_t IT_0489 = IT_0133 + IT_0303 + IT_0372 + IT_0373;
    const ccomplex_t IT_0490 = (-2)*IT_0136 + -IT_0302 + -IT_0367 + -IT_0371;
    const ccomplex_t IT_0491 = IT_0489 + IT_0490;
    const ccomplex_t IT_0492 = IT_0124*IT_0491;
    const ccomplex_t IT_0493 = IT_0178 + IT_0182 + IT_0184 + IT_0187 + IT_0188;
    const ccomplex_t IT_0494 = IT_0177*IT_0493;
    const ccomplex_t IT_0495 = IT_0342 + IT_0344 + IT_0346 + IT_0478 + IT_0480
       + IT_0484 + IT_0488 + IT_0492 + IT_0494;
    const ccomplex_t IT_0496 = IT_0169*IT_0477;
    const ccomplex_t IT_0497 = IT_0219*IT_0479;
    const ccomplex_t IT_0498 = IT_0079 + IT_0351;
    const ccomplex_t IT_0499 = -IT_0073 + -IT_0352;
    const ccomplex_t IT_0500 = IT_0498 + IT_0499;
    const ccomplex_t IT_0501 = IT_0272*IT_0500;
    const ccomplex_t IT_0502 = IT_0151 + IT_0155 + IT_0157 + IT_0160 + IT_0161
       + IT_0162;
    const ccomplex_t IT_0503 = IT_0251*IT_0502;
    const ccomplex_t IT_0504 = IT_0178 + IT_0182 + IT_0184 + IT_0187 + IT_0188
       + IT_0189;
    const ccomplex_t IT_0505 = IT_0257*IT_0504;
    const ccomplex_t IT_0506 = IT_0196*IT_0493;
    const ccomplex_t IT_0507 = IT_0203 + IT_0207 + IT_0211 + IT_0212 + IT_0213
       + IT_0214;
    const ccomplex_t IT_0508 = IT_0263*IT_0507;
    const ccomplex_t IT_0509 = IT_0099 + IT_0361;
    const ccomplex_t IT_0510 = -IT_0112 + -IT_0362;
    const ccomplex_t IT_0511 = IT_0509 + IT_0510;
    const ccomplex_t IT_0512 = IT_0281*IT_0511;
    const ccomplex_t IT_0513 = IT_0133 + IT_0371;
    const ccomplex_t IT_0514 = -IT_0136 + -IT_0373;
    const ccomplex_t IT_0515 = IT_0513 + IT_0514;
    const ccomplex_t IT_0516 = IT_0299*IT_0515;
    const ccomplex_t IT_0517 = IT_0308*IT_0500;
    const ccomplex_t IT_0518 = IT_0311*IT_0511;
    const ccomplex_t IT_0519 = IT_0317*IT_0515;
    const ccomplex_t IT_0520 = -IT_0334 + -IT_0337 + -IT_0340 + -IT_0341 + 
      -IT_0343 + -IT_0345 + 0.5*IT_0398 + 0.5*IT_0402 + 0.5*IT_0406 + (-0.5)
      *IT_0410 + (-0.5)*IT_0412 + (-0.5)*IT_0414 + -IT_0496 + -IT_0497 + (-0.5)
      *IT_0501 + 2*IT_0503 + 2*IT_0505 + -IT_0506 + 2*IT_0508 + (-0.5)*IT_0512 +
       (-0.5)*IT_0516 + 0.5*IT_0517 + 0.5*IT_0518 + 0.5*IT_0519;
    const ccomplex_t IT_0521 = IT_0495 + IT_0520;
    const ccomplex_t IT_0522 = IT_0005*IT_0521;
    const ccomplex_t IT_0523 = 1.4142135623731*IT_0522;
    const ccomplex_t IT_0524 = 0.25*IT_0523;
    return create_ccomplex_return((0 + _Complex_I*0.25)*IT_0326 + (0 +
       _Complex_I*0.25)*IT_0419 + (0 + _Complex_I*0.25)*IT_0476 + (0 +
       _Complex_I*0.25)*IT_0524);
}

