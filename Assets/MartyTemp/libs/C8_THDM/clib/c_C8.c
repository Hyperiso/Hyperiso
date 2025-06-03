#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "marty/core/looptools_interface.h"
#include "c_C8.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_C8(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t G_F = param->G_F;
    const creal_t M_W = param->M_W;
    const creal_t g_s = param->g_s;
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
    const ccomplex_t IT_0001 = pow(g_s, -1);
    const ccomplex_t IT_0002 = pow(m_b, -1);
    const ccomplex_t IT_0003 = pow(V_tb, -1);
    const ccomplex_t IT_0004 = cpow(conj(V_ts), -1);
    const ccomplex_t IT_0005 = 9.86960440108936*IT_0000*IT_0001*IT_0002
      *IT_0003*IT_0004;
    const ccomplex_t IT_0006 = 0.101321183642338*m_c;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0008 = pow(M_W, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0008*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0008*IT_0010;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0012*IT_0014;
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = pow(m_b, 2);
    const ccomplex_t IT_0018 = pow(m_s, 2);
    const ccomplex_t IT_0019 = pow(m_c, 2);
    const ccomplex_t IT_0020 = pow(M_W, 2);
    const ccomplex_t IT_0021 = mtylt_C0iC(3, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0019, IT_0020, IT_0019, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0022 = mtylt_C0iC(6, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0019, IT_0020, IT_0019, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0023 = mtylt_C0iC(0, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0019, IT_0020, IT_0019, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0024 = IT_0021 + IT_0022 + IT_0023;
    const ccomplex_t IT_0025 = IT_0016*IT_0024;
    const ccomplex_t IT_0026 = cos(beta);
    const ccomplex_t IT_0027 = sin(beta);
    const ccomplex_t IT_0028 = cpow(IT_0027, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0007*IT_0030*IT_0032;
    const ccomplex_t IT_0034 = 0.101321183642338*m_t;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = pow(m_t, 2);
    const ccomplex_t IT_0037 = pow(m_Hp, 2);
    const ccomplex_t IT_0038 = mtylt_C0iC(0, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0036, IT_0036, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0039 = mtylt_C0iC(3, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0036, IT_0036, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0040 = mtylt_C0iC(6, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0036, IT_0036, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0041 = IT_0038 + IT_0039 + IT_0040;
    const ccomplex_t IT_0042 = IT_0035*IT_0041;
    const ccomplex_t IT_0043 = 0.101321183642338*m_u;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*m_u*V_ub*e_em
      *IT_0008*IT_0010;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0008*IT_0010;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0007*IT_0045*IT_0047;
    const ccomplex_t IT_0049 = IT_0043*IT_0048;
    const ccomplex_t IT_0050 = pow(m_u, 2);
    const ccomplex_t IT_0051 = mtylt_C0iC(3, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0050, IT_0020, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0052 = mtylt_C0iC(6, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0050, IT_0020, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0053 = mtylt_C0iC(0, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0050, IT_0020, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0054 = IT_0051 + IT_0052 + IT_0053;
    const ccomplex_t IT_0055 = IT_0049*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_u*V_ub*e_em
      *IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = IT_0007*IT_0057*IT_0059;
    const ccomplex_t IT_0061 = IT_0043*IT_0060;
    const ccomplex_t IT_0062 = mtylt_C0iC(0, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0050, IT_0050, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0063 = mtylt_C0iC(3, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0050, IT_0050, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0064 = mtylt_C0iC(6, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0050, IT_0050, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0065 = IT_0062 + IT_0063 + IT_0064;
    const ccomplex_t IT_0066 = IT_0061*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = IT_0007*IT_0068*IT_0070;
    const ccomplex_t IT_0072 = IT_0006*IT_0071;
    const ccomplex_t IT_0073 = mtylt_C0iC(0, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0019, IT_0019, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0074 = mtylt_C0iC(3, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0019, IT_0019, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0075 = mtylt_C0iC(6, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0019, IT_0019, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0076 = IT_0073 + IT_0074 + IT_0075;
    const ccomplex_t IT_0077 = IT_0072*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0008*IT_0010;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0008*IT_0010;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = IT_0007*IT_0079*IT_0081;
    const ccomplex_t IT_0083 = IT_0034*IT_0082;
    const ccomplex_t IT_0084 = mtylt_C0iC(3, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0036, IT_0020, IT_0036, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0085 = mtylt_C0iC(6, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0036, IT_0020, IT_0036, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0086 = mtylt_C0iC(0, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0036, IT_0020, IT_0036, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0087 = IT_0084 + IT_0085 + IT_0086;
    const ccomplex_t IT_0088 = IT_0083*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = IT_0007*IT_0057*IT_0090;
    const ccomplex_t IT_0092 = 0.101321183642338*IT_0091;
    const ccomplex_t IT_0093 = m_s*IT_0063;
    const ccomplex_t IT_0094 = mtylt_C0iC(12, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0050, IT_0050, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0095 = m_s*IT_0094;
    const ccomplex_t IT_0096 = mtylt_C0iC(15, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0050, IT_0050, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0097 = m_s*IT_0096;
    const ccomplex_t IT_0098 = m_b*IT_0063;
    const ccomplex_t IT_0099 = m_b*IT_0064;
    const ccomplex_t IT_0100 = m_b*IT_0094;
    const ccomplex_t IT_0101 = mtylt_C0iC(18, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0050, IT_0050, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0102 = m_b*IT_0101;
    const ccomplex_t IT_0103 = IT_0093 + IT_0095 + IT_0097 + IT_0098 + IT_0099
       + IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = m_b*IT_0096;
    const ccomplex_t IT_0105 = 2*IT_0104;
    const ccomplex_t IT_0106 = IT_0103 + IT_0105;
    const ccomplex_t IT_0107 = IT_0092*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0008*IT_0010;
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = IT_0007*IT_0045*IT_0109;
    const ccomplex_t IT_0111 = 0.101321183642338*IT_0110;
    const ccomplex_t IT_0112 = mtylt_C0iC(15, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0050, IT_0020, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0113 = m_b*IT_0112;
    const ccomplex_t IT_0114 = 2*IT_0113;
    const ccomplex_t IT_0115 = m_s*IT_0052;
    const ccomplex_t IT_0116 = m_s*IT_0112;
    const ccomplex_t IT_0117 = mtylt_C0iC(18, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0050, IT_0020, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0118 = m_s*IT_0117;
    const ccomplex_t IT_0119 = m_b*IT_0051;
    const ccomplex_t IT_0120 = m_b*IT_0052;
    const ccomplex_t IT_0121 = mtylt_C0iC(12, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0050, IT_0020, IT_0050, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0122 = m_b*IT_0121;
    const ccomplex_t IT_0123 = m_b*IT_0117;
    const ccomplex_t IT_0124 = IT_0115 + IT_0116 + IT_0118 + IT_0119 + IT_0120
       + IT_0122 + IT_0123;
    const ccomplex_t IT_0125 = IT_0114 + IT_0124;
    const ccomplex_t IT_0126 = IT_0111*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = IT_0007*IT_0068*IT_0128;
    const ccomplex_t IT_0130 = 0.101321183642338*IT_0129;
    const ccomplex_t IT_0131 = m_s*IT_0074;
    const ccomplex_t IT_0132 = mtylt_C0iC(12, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0019, IT_0019, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0133 = m_s*IT_0132;
    const ccomplex_t IT_0134 = mtylt_C0iC(15, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0019, IT_0019, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0135 = m_s*IT_0134;
    const ccomplex_t IT_0136 = m_b*IT_0074;
    const ccomplex_t IT_0137 = m_b*IT_0075;
    const ccomplex_t IT_0138 = m_b*IT_0132;
    const ccomplex_t IT_0139 = mtylt_C0iC(18, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0019, IT_0019, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0140 = m_b*IT_0139;
    const ccomplex_t IT_0141 = IT_0131 + IT_0133 + IT_0135 + IT_0136 + IT_0137
       + IT_0138 + IT_0140;
    const ccomplex_t IT_0142 = m_b*IT_0134;
    const ccomplex_t IT_0143 = 2*IT_0142;
    const ccomplex_t IT_0144 = IT_0141 + IT_0143;
    const ccomplex_t IT_0145 = IT_0130*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0008*IT_0010;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = IT_0007*IT_0012*IT_0147;
    const ccomplex_t IT_0149 = 0.101321183642338*IT_0148;
    const ccomplex_t IT_0150 = mtylt_C0iC(15, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0019, IT_0020, IT_0019, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0151 = m_b*IT_0150;
    const ccomplex_t IT_0152 = 2*IT_0151;
    const ccomplex_t IT_0153 = m_s*IT_0022;
    const ccomplex_t IT_0154 = m_s*IT_0150;
    const ccomplex_t IT_0155 = mtylt_C0iC(18, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0019, IT_0020, IT_0019, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0156 = m_s*IT_0155;
    const ccomplex_t IT_0157 = m_b*IT_0021;
    const ccomplex_t IT_0158 = m_b*IT_0022;
    const ccomplex_t IT_0159 = mtylt_C0iC(12, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0019, IT_0020, IT_0019, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0160 = m_b*IT_0159;
    const ccomplex_t IT_0161 = m_b*IT_0155;
    const ccomplex_t IT_0162 = IT_0153 + IT_0154 + IT_0156 + IT_0157 + IT_0158
       + IT_0160 + IT_0161;
    const ccomplex_t IT_0163 = IT_0152 + IT_0162;
    const ccomplex_t IT_0164 = IT_0149*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0166 = 0.5*IT_0165;
    const ccomplex_t IT_0167 = IT_0007*IT_0030*IT_0166;
    const ccomplex_t IT_0168 = 0.101321183642338*IT_0167;
    const ccomplex_t IT_0169 = m_s*IT_0039;
    const ccomplex_t IT_0170 = mtylt_C0iC(12, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0036, IT_0036, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0171 = m_s*IT_0170;
    const ccomplex_t IT_0172 = mtylt_C0iC(15, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0036, IT_0036, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0173 = m_s*IT_0172;
    const ccomplex_t IT_0174 = m_b*IT_0039;
    const ccomplex_t IT_0175 = m_b*IT_0040;
    const ccomplex_t IT_0176 = m_b*IT_0170;
    const ccomplex_t IT_0177 = mtylt_C0iC(18, IT_0017 + IT_0018 + (-2)*s_12,
       IT_0018, IT_0017, IT_0036, IT_0036, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0178 = m_b*IT_0177;
    const ccomplex_t IT_0179 = IT_0169 + IT_0171 + IT_0173 + IT_0174 + IT_0175
       + IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = m_b*IT_0172;
    const ccomplex_t IT_0181 = 2*IT_0180;
    const ccomplex_t IT_0182 = IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = IT_0168*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0008*IT_0010;
    const ccomplex_t IT_0185 = 0.5*IT_0184;
    const ccomplex_t IT_0186 = IT_0007*IT_0079*IT_0185;
    const ccomplex_t IT_0187 = 0.101321183642338*IT_0186;
    const ccomplex_t IT_0188 = mtylt_C0iC(15, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0036, IT_0020, IT_0036, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0189 = m_b*IT_0188;
    const ccomplex_t IT_0190 = 2*IT_0189;
    const ccomplex_t IT_0191 = m_s*IT_0085;
    const ccomplex_t IT_0192 = m_s*IT_0188;
    const ccomplex_t IT_0193 = mtylt_C0iC(18, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0036, IT_0020, IT_0036, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0194 = m_s*IT_0193;
    const ccomplex_t IT_0195 = m_b*IT_0084;
    const ccomplex_t IT_0196 = m_b*IT_0085;
    const ccomplex_t IT_0197 = mtylt_C0iC(12, IT_0017, IT_0018, IT_0017 +
       IT_0018 + (-2)*s_12, IT_0036, IT_0020, IT_0036, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0198 = m_b*IT_0197;
    const ccomplex_t IT_0199 = m_b*IT_0193;
    const ccomplex_t IT_0200 = IT_0191 + IT_0192 + IT_0194 + IT_0195 + IT_0196
       + IT_0198 + IT_0199;
    const ccomplex_t IT_0201 = IT_0190 + IT_0200;
    const ccomplex_t IT_0202 = IT_0187*IT_0201;
    const ccomplex_t IT_0203 = IT_0025 + IT_0042 + IT_0055 + IT_0066 + IT_0077
       + IT_0088 + IT_0107 + IT_0126 + IT_0145 + IT_0164 + IT_0183 + IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0205 = (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = IT_0007*IT_0128*IT_0205;
    const ccomplex_t IT_0207 = IT_0006*IT_0206;
    const ccomplex_t IT_0208 = IT_0076*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0008*IT_0010;
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0007*IT_0185*IT_0210;
    const ccomplex_t IT_0212 = IT_0034*IT_0211;
    const ccomplex_t IT_0213 = IT_0087*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*m_b*V_ub*e_em
      *IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0215 = (-0.5)*IT_0214;
    const ccomplex_t IT_0216 = IT_0007*IT_0090*IT_0215;
    const ccomplex_t IT_0217 = IT_0043*IT_0216;
    const ccomplex_t IT_0218 = IT_0065*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1.4142135623731)*m_b*V_ub*e_em
      *IT_0008*IT_0010;
    const ccomplex_t IT_0220 = (-0.5)*IT_0219;
    const ccomplex_t IT_0221 = IT_0007*IT_0109*IT_0220;
    const ccomplex_t IT_0222 = IT_0043*IT_0221;
    const ccomplex_t IT_0223 = IT_0054*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0008*IT_0010;
    const ccomplex_t IT_0225 = (-0.5)*IT_0224;
    const ccomplex_t IT_0226 = IT_0007*IT_0147*IT_0225;
    const ccomplex_t IT_0227 = IT_0006*IT_0226;
    const ccomplex_t IT_0228 = IT_0024*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0008*IT_0010*IT_0026*IT_0028;
    const ccomplex_t IT_0230 = (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = IT_0007*IT_0166*IT_0230;
    const ccomplex_t IT_0232 = IT_0034*IT_0231;
    const ccomplex_t IT_0233 = IT_0041*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0010;
    const ccomplex_t IT_0235 = 0.5*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1.4142135623731)*V_ub*e_em
      *IT_0010;
    const ccomplex_t IT_0237 = 0.5*IT_0236;
    const ccomplex_t IT_0238 = IT_0007*IT_0235*IT_0237;
    const ccomplex_t IT_0239 = 0.101321183642338*IT_0238;
    const ccomplex_t IT_0240 = m_s*IT_0051;
    const ccomplex_t IT_0241 = IT_0115 + IT_0116 + IT_0118 + IT_0119 + IT_0120
       + IT_0122 + IT_0123 + IT_0240;
    const ccomplex_t IT_0242 = IT_0114 + IT_0241;
    const ccomplex_t IT_0243 = IT_0239*IT_0242;
    const ccomplex_t IT_0244 = IT_0007*IT_0059*IT_0215;
    const ccomplex_t IT_0245 = 0.101321183642338*IT_0244;
    const ccomplex_t IT_0246 = IT_0106*IT_0245;
    const ccomplex_t IT_0247 = IT_0007*IT_0047*IT_0220;
    const ccomplex_t IT_0248 = 0.101321183642338*IT_0247;
    const ccomplex_t IT_0249 = IT_0125*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0010;
    const ccomplex_t IT_0251 = 0.5*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0010;
    const ccomplex_t IT_0253 = 0.5*IT_0252;
    const ccomplex_t IT_0254 = IT_0007*IT_0251*IT_0253;
    const ccomplex_t IT_0255 = 0.101321183642338*IT_0254;
    const ccomplex_t IT_0256 = m_s*IT_0021;
    const ccomplex_t IT_0257 = IT_0153 + IT_0154 + IT_0156 + IT_0157 + IT_0158
       + IT_0160 + IT_0161 + IT_0256;
    const ccomplex_t IT_0258 = IT_0152 + IT_0257;
    const ccomplex_t IT_0259 = IT_0255*IT_0258;
    const ccomplex_t IT_0260 = IT_0007*IT_0070*IT_0205;
    const ccomplex_t IT_0261 = 0.101321183642338*IT_0260;
    const ccomplex_t IT_0262 = IT_0144*IT_0261;
    const ccomplex_t IT_0263 = IT_0007*IT_0014*IT_0225;
    const ccomplex_t IT_0264 = 0.101321183642338*IT_0263;
    const ccomplex_t IT_0265 = IT_0163*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0010;
    const ccomplex_t IT_0267 = 0.5*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0010;
    const ccomplex_t IT_0269 = 0.5*IT_0268;
    const ccomplex_t IT_0270 = IT_0007*IT_0267*IT_0269;
    const ccomplex_t IT_0271 = 0.101321183642338*IT_0270;
    const ccomplex_t IT_0272 = m_s*IT_0084;
    const ccomplex_t IT_0273 = IT_0191 + IT_0192 + IT_0194 + IT_0195 + IT_0196
       + IT_0198 + IT_0199 + IT_0272;
    const ccomplex_t IT_0274 = IT_0190 + IT_0273;
    const ccomplex_t IT_0275 = IT_0271*IT_0274;
    const ccomplex_t IT_0276 = IT_0007*IT_0032*IT_0230;
    const ccomplex_t IT_0277 = 0.101321183642338*IT_0276;
    const ccomplex_t IT_0278 = IT_0182*IT_0277;
    const ccomplex_t IT_0279 = IT_0007*IT_0081*IT_0210;
    const ccomplex_t IT_0280 = 0.101321183642338*IT_0279;
    const ccomplex_t IT_0281 = IT_0201*IT_0280;
    const ccomplex_t IT_0282 = -IT_0208 + -IT_0213 + -IT_0218 + -IT_0223 + 
      -IT_0228 + -IT_0233 + 2*IT_0243 + -IT_0246 + -IT_0249 + 2*IT_0259 + 
      -IT_0262 + -IT_0265 + 2*IT_0275 + -IT_0278 + -IT_0281;
    const ccomplex_t IT_0283 = IT_0203 + IT_0282;
    const ccomplex_t IT_0284 = IT_0005*IT_0283;
    const ccomplex_t IT_0285 = 1.4142135623731*IT_0284;
    const ccomplex_t IT_0286 = (-0.25)*IT_0285;
    const ccomplex_t IT_0287 = IT_0153 + IT_0154 + IT_0156 + IT_0256;
    const ccomplex_t IT_0288 = (-2)*IT_0151 + -IT_0157 + -IT_0158 + -IT_0160 +
       -IT_0161;
    const ccomplex_t IT_0289 = IT_0287 + IT_0288;
    const ccomplex_t IT_0290 = IT_0255*IT_0289;
    const ccomplex_t IT_0291 = IT_0191 + IT_0192 + IT_0194 + IT_0272;
    const ccomplex_t IT_0292 = (-2)*IT_0189 + -IT_0195 + -IT_0196 + -IT_0198 +
       -IT_0199;
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = IT_0271*IT_0293;
    const ccomplex_t IT_0295 = IT_0115 + IT_0116 + IT_0118 + IT_0240;
    const ccomplex_t IT_0296 = (-2)*IT_0113 + -IT_0119 + -IT_0120 + -IT_0122 +
       -IT_0123;
    const ccomplex_t IT_0297 = IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = IT_0239*IT_0297;
    const ccomplex_t IT_0299 = IT_0290 + IT_0294 + IT_0298;
    const ccomplex_t IT_0300 = IT_0093 + IT_0095 + IT_0097;
    const ccomplex_t IT_0301 = -IT_0098 + -IT_0099 + -IT_0100 + -IT_0102 + (-2
      )*IT_0104;
    const ccomplex_t IT_0302 = IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = IT_0092*IT_0302;
    const ccomplex_t IT_0304 = IT_0115 + IT_0116 + IT_0118;
    const ccomplex_t IT_0305 = IT_0296 + IT_0304;
    const ccomplex_t IT_0306 = IT_0111*IT_0305;
    const ccomplex_t IT_0307 = IT_0245*IT_0302;
    const ccomplex_t IT_0308 = IT_0248*IT_0305;
    const ccomplex_t IT_0309 = IT_0131 + IT_0133 + IT_0135;
    const ccomplex_t IT_0310 = -IT_0136 + -IT_0137 + -IT_0138 + -IT_0140 + (-2
      )*IT_0142;
    const ccomplex_t IT_0311 = IT_0309 + IT_0310;
    const ccomplex_t IT_0312 = IT_0130*IT_0311;
    const ccomplex_t IT_0313 = IT_0153 + IT_0154 + IT_0156;
    const ccomplex_t IT_0314 = IT_0288 + IT_0313;
    const ccomplex_t IT_0315 = IT_0149*IT_0314;
    const ccomplex_t IT_0316 = IT_0261*IT_0311;
    const ccomplex_t IT_0317 = IT_0264*IT_0314;
    const ccomplex_t IT_0318 = IT_0169 + IT_0171 + IT_0173;
    const ccomplex_t IT_0319 = -IT_0174 + -IT_0175 + -IT_0176 + -IT_0178 + (-2
      )*IT_0180;
    const ccomplex_t IT_0320 = IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = IT_0168*IT_0320;
    const ccomplex_t IT_0322 = IT_0191 + IT_0192 + IT_0194;
    const ccomplex_t IT_0323 = IT_0292 + IT_0322;
    const ccomplex_t IT_0324 = IT_0187*IT_0323;
    const ccomplex_t IT_0325 = IT_0277*IT_0320;
    const ccomplex_t IT_0326 = IT_0280*IT_0323;
    const ccomplex_t IT_0327 = 0.5*IT_0025 + 0.5*IT_0042 + 0.5*IT_0055 + 0.5
      *IT_0066 + 0.5*IT_0077 + 0.5*IT_0088 + 0.5*IT_0208 + 0.5*IT_0213 + 0.5
      *IT_0218 + 0.5*IT_0223 + 0.5*IT_0228 + 0.5*IT_0233 + 0.5*IT_0303 + 0.5
      *IT_0306 + 0.5*IT_0307 + 0.5*IT_0308 + 0.5*IT_0312 + 0.5*IT_0315 + 0.5
      *IT_0316 + 0.5*IT_0317 + 0.5*IT_0321 + 0.5*IT_0324 + 0.5*IT_0325 + 0.5
      *IT_0326;
    const ccomplex_t IT_0328 = IT_0299 + IT_0327;
    const ccomplex_t IT_0329 = IT_0005*IT_0328;
    const ccomplex_t IT_0330 = 1.4142135623731*IT_0329;
    const ccomplex_t IT_0331 = 0.5*IT_0330;
    const ccomplex_t IT_0332 = IT_0115 + IT_0118;
    const ccomplex_t IT_0333 = -IT_0113 + -IT_0119 + -IT_0120 + -IT_0123;
    const ccomplex_t IT_0334 = IT_0332 + IT_0333;
    const ccomplex_t IT_0335 = IT_0239*IT_0334;
    const ccomplex_t IT_0336 = IT_0153 + IT_0156;
    const ccomplex_t IT_0337 = -IT_0151 + -IT_0157 + -IT_0158 + -IT_0161;
    const ccomplex_t IT_0338 = IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = IT_0255*IT_0338;
    const ccomplex_t IT_0340 = IT_0191 + IT_0194;
    const ccomplex_t IT_0341 = -IT_0189 + -IT_0195 + -IT_0196 + -IT_0199;
    const ccomplex_t IT_0342 = IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = IT_0271*IT_0342;
    const ccomplex_t IT_0344 = IT_0335 + IT_0339 + IT_0343;
    const ccomplex_t IT_0345 = IT_0061*IT_0063;
    const ccomplex_t IT_0346 = IT_0049*IT_0052;
    const ccomplex_t IT_0347 = IT_0063*IT_0217;
    const ccomplex_t IT_0348 = IT_0052*IT_0222;
    const ccomplex_t IT_0349 = IT_0072*IT_0074;
    const ccomplex_t IT_0350 = IT_0016*IT_0022;
    const ccomplex_t IT_0351 = IT_0074*IT_0207;
    const ccomplex_t IT_0352 = IT_0022*IT_0227;
    const ccomplex_t IT_0353 = IT_0035*IT_0039;
    const ccomplex_t IT_0354 = IT_0083*IT_0085;
    const ccomplex_t IT_0355 = IT_0039*IT_0232;
    const ccomplex_t IT_0356 = IT_0085*IT_0212;
    const ccomplex_t IT_0357 = IT_0093 + IT_0095;
    const ccomplex_t IT_0358 = -IT_0098 + -IT_0100 + -IT_0104;
    const ccomplex_t IT_0359 = IT_0357 + IT_0358;
    const ccomplex_t IT_0360 = IT_0092*IT_0359;
    const ccomplex_t IT_0361 = -IT_0113 + -IT_0120 + -IT_0123;
    const ccomplex_t IT_0362 = IT_0332 + IT_0361;
    const ccomplex_t IT_0363 = IT_0111*IT_0362;
    const ccomplex_t IT_0364 = IT_0245*IT_0359;
    const ccomplex_t IT_0365 = IT_0248*IT_0362;
    const ccomplex_t IT_0366 = IT_0131 + IT_0133;
    const ccomplex_t IT_0367 = -IT_0136 + -IT_0138 + -IT_0142;
    const ccomplex_t IT_0368 = IT_0366 + IT_0367;
    const ccomplex_t IT_0369 = IT_0130*IT_0368;
    const ccomplex_t IT_0370 = -IT_0151 + -IT_0158 + -IT_0161;
    const ccomplex_t IT_0371 = IT_0336 + IT_0370;
    const ccomplex_t IT_0372 = IT_0149*IT_0371;
    const ccomplex_t IT_0373 = IT_0261*IT_0368;
    const ccomplex_t IT_0374 = IT_0264*IT_0371;
    const ccomplex_t IT_0375 = -IT_0189 + -IT_0196 + -IT_0199;
    const ccomplex_t IT_0376 = IT_0340 + IT_0375;
    const ccomplex_t IT_0377 = IT_0280*IT_0376;
    const ccomplex_t IT_0378 = IT_0169 + IT_0171;
    const ccomplex_t IT_0379 = -IT_0174 + -IT_0176 + -IT_0180;
    const ccomplex_t IT_0380 = IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = IT_0168*IT_0380;
    const ccomplex_t IT_0382 = IT_0187*IT_0376;
    const ccomplex_t IT_0383 = IT_0277*IT_0380;
    const ccomplex_t IT_0384 = 0.5*IT_0345 + 0.5*IT_0346 + 0.5*IT_0347 + 0.5
      *IT_0348 + 0.5*IT_0349 + 0.5*IT_0350 + 0.5*IT_0351 + 0.5*IT_0352 + 0.5
      *IT_0353 + 0.5*IT_0354 + 0.5*IT_0355 + 0.5*IT_0356 + 0.5*IT_0360 + 0.5
      *IT_0363 + 0.5*IT_0364 + 0.5*IT_0365 + 0.5*IT_0369 + 0.5*IT_0372 + 0.5
      *IT_0373 + 0.5*IT_0374 + 0.5*IT_0377 + 0.5*IT_0381 + 0.5*IT_0382 + 0.5
      *IT_0383;
    const ccomplex_t IT_0385 = IT_0344 + IT_0384;
    const ccomplex_t IT_0386 = IT_0005*IT_0385;
    const ccomplex_t IT_0387 = 1.4142135623731*IT_0386;
    const ccomplex_t IT_0388 = (-0.5)*IT_0387;
    const ccomplex_t IT_0389 = IT_0093 + IT_0095 + IT_0098 + IT_0100 + IT_0104;
    const ccomplex_t IT_0390 = IT_0092*IT_0389;
    const ccomplex_t IT_0391 = IT_0113 + IT_0115 + IT_0118 + IT_0120 + IT_0123;
    const ccomplex_t IT_0392 = IT_0111*IT_0391;
    const ccomplex_t IT_0393 = IT_0131 + IT_0133 + IT_0136 + IT_0138 + IT_0142;
    const ccomplex_t IT_0394 = IT_0130*IT_0393;
    const ccomplex_t IT_0395 = IT_0151 + IT_0153 + IT_0156 + IT_0158 + IT_0161;
    const ccomplex_t IT_0396 = IT_0149*IT_0395;
    const ccomplex_t IT_0397 = IT_0169 + IT_0171 + IT_0174 + IT_0176 + IT_0180;
    const ccomplex_t IT_0398 = IT_0168*IT_0397;
    const ccomplex_t IT_0399 = IT_0189 + IT_0191 + IT_0194 + IT_0196 + IT_0199;
    const ccomplex_t IT_0400 = IT_0187*IT_0399;
    const ccomplex_t IT_0401 = IT_0345 + IT_0346 + IT_0349 + IT_0350 + IT_0353
       + IT_0354 + IT_0390 + IT_0392 + IT_0394 + IT_0396 + IT_0398 + IT_0400;
    const ccomplex_t IT_0402 = IT_0113 + IT_0115 + IT_0118 + IT_0119 + IT_0120
       + IT_0123;
    const ccomplex_t IT_0403 = IT_0239*IT_0402;
    const ccomplex_t IT_0404 = IT_0245*IT_0389;
    const ccomplex_t IT_0405 = IT_0248*IT_0391;
    const ccomplex_t IT_0406 = IT_0151 + IT_0153 + IT_0156 + IT_0157 + IT_0158
       + IT_0161;
    const ccomplex_t IT_0407 = IT_0255*IT_0406;
    const ccomplex_t IT_0408 = IT_0261*IT_0393;
    const ccomplex_t IT_0409 = IT_0264*IT_0395;
    const ccomplex_t IT_0410 = IT_0189 + IT_0191 + IT_0194 + IT_0195 + IT_0196
       + IT_0199;
    const ccomplex_t IT_0411 = IT_0271*IT_0410;
    const ccomplex_t IT_0412 = IT_0277*IT_0397;
    const ccomplex_t IT_0413 = IT_0280*IT_0399;
    const ccomplex_t IT_0414 = -IT_0347 + -IT_0348 + -IT_0351 + -IT_0352 + 
      -IT_0355 + -IT_0356 + 2*IT_0403 + -IT_0404 + -IT_0405 + 2*IT_0407 + 
      -IT_0408 + -IT_0409 + 2*IT_0411 + -IT_0412 + -IT_0413;
    const ccomplex_t IT_0415 = IT_0401 + IT_0414;
    const ccomplex_t IT_0416 = IT_0005*IT_0415;
    const ccomplex_t IT_0417 = 1.4142135623731*IT_0416;
    const ccomplex_t IT_0418 = 0.25*IT_0417;
    return create_ccomplex_return((0 + _Complex_I*0.25)*IT_0286 + (0 +
       _Complex_I*0.25)*IT_0331 + (0 + _Complex_I*0.25)*IT_0388 + (0 +
       _Complex_I*0.25)*IT_0418);
}

