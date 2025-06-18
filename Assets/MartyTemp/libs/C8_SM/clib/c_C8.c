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
    const creal_t e_em = param->e_em;
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
    const ccomplex_t IT_0006 = pow(m_b, 2);
    const ccomplex_t IT_0007 = pow(m_s, 2);
    const ccomplex_t IT_0008 = pow(m_u, 2);
    const ccomplex_t IT_0009 = pow(M_W, 2);
    const ccomplex_t IT_0010 = mtylt_C0iC(3, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0011 = m_s*IT_0010;
    const ccomplex_t IT_0012 = mtylt_C0iC(6, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0013 = m_s*IT_0012;
    const ccomplex_t IT_0014 = mtylt_C0iC(15, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0015 = m_s*IT_0014;
    const ccomplex_t IT_0016 = mtylt_C0iC(18, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0017 = m_s*IT_0016;
    const ccomplex_t IT_0018 = IT_0011 + IT_0013 + IT_0015 + IT_0017;
    const ccomplex_t IT_0019 = m_b*IT_0010;
    const ccomplex_t IT_0020 = m_b*IT_0012;
    const ccomplex_t IT_0021 = mtylt_C0iC(12, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0022 = m_b*IT_0021;
    const ccomplex_t IT_0023 = m_b*IT_0014;
    const ccomplex_t IT_0024 = m_b*IT_0016;
    const ccomplex_t IT_0025 = -IT_0019 + -IT_0020 + -IT_0022 + (-2)*IT_0023 +
       -IT_0024;
    const ccomplex_t IT_0026 = IT_0018 + IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0028 = sin(theta_W);
    const ccomplex_t IT_0029 = cpow(IT_0028, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*V_ub*e_em
      *IT_0029;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0031*IT_0033;
    const ccomplex_t IT_0035 = 0.101321183642338*IT_0034;
    const ccomplex_t IT_0036 = IT_0026*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0029;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0029;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0027*IT_0038*IT_0040;
    const ccomplex_t IT_0042 = 0.101321183642338*IT_0041;
    const ccomplex_t IT_0043 = pow(m_c, 2);
    const ccomplex_t IT_0044 = mtylt_C0iC(3, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0043, IT_0009, IT_0043, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0045 = m_s*IT_0044;
    const ccomplex_t IT_0046 = mtylt_C0iC(6, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0043, IT_0009, IT_0043, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0047 = m_s*IT_0046;
    const ccomplex_t IT_0048 = mtylt_C0iC(15, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0043, IT_0009, IT_0043, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0049 = m_s*IT_0048;
    const ccomplex_t IT_0050 = mtylt_C0iC(18, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0043, IT_0009, IT_0043, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0051 = m_s*IT_0050;
    const ccomplex_t IT_0052 = IT_0045 + IT_0047 + IT_0049 + IT_0051;
    const ccomplex_t IT_0053 = m_b*IT_0044;
    const ccomplex_t IT_0054 = m_b*IT_0046;
    const ccomplex_t IT_0055 = mtylt_C0iC(12, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0043, IT_0009, IT_0043, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0056 = m_b*IT_0055;
    const ccomplex_t IT_0057 = m_b*IT_0048;
    const ccomplex_t IT_0058 = m_b*IT_0050;
    const ccomplex_t IT_0059 = -IT_0053 + -IT_0054 + -IT_0056 + (-2)*IT_0057 +
       -IT_0058;
    const ccomplex_t IT_0060 = IT_0052 + IT_0059;
    const ccomplex_t IT_0061 = IT_0042*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0029;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0029;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0027*IT_0063*IT_0065;
    const ccomplex_t IT_0067 = 0.101321183642338*IT_0066;
    const ccomplex_t IT_0068 = pow(m_t, 2);
    const ccomplex_t IT_0069 = mtylt_C0iC(3, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0068, IT_0009, IT_0068, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0070 = m_s*IT_0069;
    const ccomplex_t IT_0071 = mtylt_C0iC(6, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0068, IT_0009, IT_0068, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0072 = m_s*IT_0071;
    const ccomplex_t IT_0073 = mtylt_C0iC(15, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0068, IT_0009, IT_0068, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0074 = m_s*IT_0073;
    const ccomplex_t IT_0075 = mtylt_C0iC(18, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0068, IT_0009, IT_0068, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0076 = m_s*IT_0075;
    const ccomplex_t IT_0077 = IT_0070 + IT_0072 + IT_0074 + IT_0076;
    const ccomplex_t IT_0078 = mtylt_C0iC(12, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0068, IT_0009, IT_0068, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0079 = m_b*IT_0078;
    const ccomplex_t IT_0080 = m_b*IT_0069;
    const ccomplex_t IT_0081 = m_b*IT_0071;
    const ccomplex_t IT_0082 = m_b*IT_0073;
    const ccomplex_t IT_0083 = m_b*IT_0075;
    const ccomplex_t IT_0084 = -IT_0079 + -IT_0080 + -IT_0081 + (-2)*IT_0082 +
       -IT_0083;
    const ccomplex_t IT_0085 = IT_0077 + IT_0084;
    const ccomplex_t IT_0086 = IT_0067*IT_0085;
    const ccomplex_t IT_0087 = IT_0036 + IT_0061 + IT_0086;
    const ccomplex_t IT_0088 = IT_0013 + IT_0015 + IT_0017;
    const ccomplex_t IT_0089 = IT_0025 + IT_0088;
    const ccomplex_t IT_0090 = pow(M_W, -1);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0029*IT_0090;
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_b*V_ub*e_em
      *IT_0029*IT_0090;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = IT_0027*IT_0092*IT_0094;
    const ccomplex_t IT_0096 = 0.101321183642338*IT_0095;
    const ccomplex_t IT_0097 = IT_0089*IT_0096;
    const ccomplex_t IT_0098 = mtylt_C0iC(0, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0068, IT_0009, IT_0068, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0099 = IT_0069 + IT_0071 + IT_0098;
    const ccomplex_t IT_0100 = 0.101321183642338*m_t;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0029*IT_0090;
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0029*IT_0090;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0027*IT_0102*IT_0104;
    const ccomplex_t IT_0106 = IT_0100*IT_0105;
    const ccomplex_t IT_0107 = IT_0099*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0029*IT_0090;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0029*IT_0090;
    const ccomplex_t IT_0111 = 0.5*IT_0110;
    const ccomplex_t IT_0112 = IT_0027*IT_0109*IT_0111;
    const ccomplex_t IT_0113 = IT_0100*IT_0112;
    const ccomplex_t IT_0114 = IT_0099*IT_0113;
    const ccomplex_t IT_0115 = mtylt_C0iC(0, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0009, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0116 = IT_0010 + IT_0012 + IT_0115;
    const ccomplex_t IT_0117 = 0.101321183642338*m_u;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*m_u*V_ub*e_em
      *IT_0029*IT_0090;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0027*IT_0092*IT_0119;
    const ccomplex_t IT_0121 = IT_0117*IT_0120;
    const ccomplex_t IT_0122 = IT_0116*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0029*IT_0090;
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = IT_0027*IT_0094*IT_0124;
    const ccomplex_t IT_0126 = IT_0117*IT_0125;
    const ccomplex_t IT_0127 = IT_0116*IT_0126;
    const ccomplex_t IT_0128 = mtylt_C0iC(0, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0043, IT_0009, IT_0043, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0129 = IT_0044 + IT_0046 + IT_0128;
    const ccomplex_t IT_0130 = 0.101321183642338*m_c;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0029*IT_0090;
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0029*IT_0090;
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = IT_0027*IT_0132*IT_0134;
    const ccomplex_t IT_0136 = IT_0130*IT_0135;
    const ccomplex_t IT_0137 = IT_0129*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0029*IT_0090;
    const ccomplex_t IT_0139 = (-0.5)*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0029*IT_0090;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0027*IT_0139*IT_0141;
    const ccomplex_t IT_0143 = IT_0130*IT_0142;
    const ccomplex_t IT_0144 = IT_0129*IT_0143;
    const ccomplex_t IT_0145 = IT_0027*IT_0119*IT_0124;
    const ccomplex_t IT_0146 = 0.101321183642338*IT_0145;
    const ccomplex_t IT_0147 = IT_0089*IT_0146;
    const ccomplex_t IT_0148 = IT_0027*IT_0134*IT_0141;
    const ccomplex_t IT_0149 = 0.101321183642338*IT_0148;
    const ccomplex_t IT_0150 = IT_0047 + IT_0049 + IT_0051;
    const ccomplex_t IT_0151 = IT_0059 + IT_0150;
    const ccomplex_t IT_0152 = IT_0149*IT_0151;
    const ccomplex_t IT_0153 = IT_0027*IT_0132*IT_0139;
    const ccomplex_t IT_0154 = 0.101321183642338*IT_0153;
    const ccomplex_t IT_0155 = IT_0151*IT_0154;
    const ccomplex_t IT_0156 = IT_0027*IT_0104*IT_0111;
    const ccomplex_t IT_0157 = 0.101321183642338*IT_0156;
    const ccomplex_t IT_0158 = IT_0072 + IT_0074 + IT_0076;
    const ccomplex_t IT_0159 = IT_0084 + IT_0158;
    const ccomplex_t IT_0160 = IT_0157*IT_0159;
    const ccomplex_t IT_0161 = IT_0027*IT_0102*IT_0109;
    const ccomplex_t IT_0162 = 0.101321183642338*IT_0161;
    const ccomplex_t IT_0163 = IT_0159*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0097 + 0.5*IT_0107 + 0.5*IT_0114 + 0.5
      *IT_0122 + 0.5*IT_0127 + 0.5*IT_0137 + 0.5*IT_0144 + 0.5*IT_0147 + 0.5
      *IT_0152 + 0.5*IT_0155 + 0.5*IT_0160 + 0.5*IT_0163;
    const ccomplex_t IT_0165 = IT_0087 + IT_0164;
    const ccomplex_t IT_0166 = IT_0005*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = 0.5*IT_0167;
    const ccomplex_t IT_0169 = IT_0013 + IT_0017;
    const ccomplex_t IT_0170 = -IT_0019 + -IT_0020 + -IT_0023 + -IT_0024;
    const ccomplex_t IT_0171 = IT_0169 + IT_0170;
    const ccomplex_t IT_0172 = IT_0035*IT_0171;
    const ccomplex_t IT_0173 = IT_0047 + IT_0051;
    const ccomplex_t IT_0174 = -IT_0053 + -IT_0054 + -IT_0057 + -IT_0058;
    const ccomplex_t IT_0175 = IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = IT_0042*IT_0175;
    const ccomplex_t IT_0177 = IT_0072 + IT_0076;
    const ccomplex_t IT_0178 = -IT_0080 + -IT_0081 + -IT_0082 + -IT_0083;
    const ccomplex_t IT_0179 = IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = IT_0067*IT_0179;
    const ccomplex_t IT_0181 = IT_0172 + IT_0176 + IT_0180;
    const ccomplex_t IT_0182 = -IT_0020 + -IT_0023 + -IT_0024;
    const ccomplex_t IT_0183 = IT_0169 + IT_0182;
    const ccomplex_t IT_0184 = IT_0096*IT_0183;
    const ccomplex_t IT_0185 = IT_0012*IT_0126;
    const ccomplex_t IT_0186 = IT_0046*IT_0136;
    const ccomplex_t IT_0187 = IT_0046*IT_0143;
    const ccomplex_t IT_0188 = IT_0071*IT_0106;
    const ccomplex_t IT_0189 = IT_0071*IT_0113;
    const ccomplex_t IT_0190 = IT_0012*IT_0121;
    const ccomplex_t IT_0191 = IT_0146*IT_0183;
    const ccomplex_t IT_0192 = -IT_0054 + -IT_0057 + -IT_0058;
    const ccomplex_t IT_0193 = IT_0173 + IT_0192;
    const ccomplex_t IT_0194 = IT_0149*IT_0193;
    const ccomplex_t IT_0195 = IT_0154*IT_0193;
    const ccomplex_t IT_0196 = -IT_0081 + -IT_0082 + -IT_0083;
    const ccomplex_t IT_0197 = IT_0177 + IT_0196;
    const ccomplex_t IT_0198 = IT_0157*IT_0197;
    const ccomplex_t IT_0199 = IT_0162*IT_0197;
    const ccomplex_t IT_0200 = 0.5*IT_0184 + 0.5*IT_0185 + 0.5*IT_0186 + 0.5
      *IT_0187 + 0.5*IT_0188 + 0.5*IT_0189 + 0.5*IT_0190 + 0.5*IT_0191 + 0.5
      *IT_0194 + 0.5*IT_0195 + 0.5*IT_0198 + 0.5*IT_0199;
    const ccomplex_t IT_0201 = IT_0181 + IT_0200;
    const ccomplex_t IT_0202 = IT_0005*IT_0201;
    const ccomplex_t IT_0203 = 1.4142135623731*IT_0202;
    const ccomplex_t IT_0204 = (-0.5)*IT_0203;
    const ccomplex_t IT_0205 = IT_0011 + IT_0013 + IT_0015 + IT_0017 + IT_0019
       + IT_0020 + IT_0022 + IT_0024;
    const ccomplex_t IT_0206 = 2*IT_0023;
    const ccomplex_t IT_0207 = IT_0205 + IT_0206;
    const ccomplex_t IT_0208 = IT_0035*IT_0207;
    const ccomplex_t IT_0209 = IT_0045 + IT_0047 + IT_0049 + IT_0051 + IT_0053
       + IT_0054 + IT_0056 + IT_0058;
    const ccomplex_t IT_0210 = 2*IT_0057;
    const ccomplex_t IT_0211 = IT_0209 + IT_0210;
    const ccomplex_t IT_0212 = IT_0042*IT_0211;
    const ccomplex_t IT_0213 = 2*IT_0082;
    const ccomplex_t IT_0214 = IT_0070 + IT_0072 + IT_0074 + IT_0076 + IT_0079
       + IT_0080 + IT_0081 + IT_0083;
    const ccomplex_t IT_0215 = IT_0213 + IT_0214;
    const ccomplex_t IT_0216 = IT_0067*IT_0215;
    const ccomplex_t IT_0217 = IT_0208 + IT_0212 + IT_0216;
    const ccomplex_t IT_0218 = IT_0013 + IT_0015 + IT_0017 + IT_0019 + IT_0020
       + IT_0022 + IT_0024;
    const ccomplex_t IT_0219 = IT_0206 + IT_0218;
    const ccomplex_t IT_0220 = IT_0146*IT_0219;
    const ccomplex_t IT_0221 = IT_0096*IT_0219;
    const ccomplex_t IT_0222 = IT_0047 + IT_0049 + IT_0051 + IT_0053 + IT_0054
       + IT_0056 + IT_0058;
    const ccomplex_t IT_0223 = IT_0210 + IT_0222;
    const ccomplex_t IT_0224 = IT_0149*IT_0223;
    const ccomplex_t IT_0225 = IT_0154*IT_0223;
    const ccomplex_t IT_0226 = IT_0072 + IT_0074 + IT_0076 + IT_0079 + IT_0080
       + IT_0081 + IT_0083;
    const ccomplex_t IT_0227 = IT_0213 + IT_0226;
    const ccomplex_t IT_0228 = IT_0157*IT_0227;
    const ccomplex_t IT_0229 = IT_0162*IT_0227;
    const ccomplex_t IT_0230 = 0.5*IT_0107 + (-0.5)*IT_0114 + 0.5*IT_0122 + (
      -0.5)*IT_0127 + 0.5*IT_0137 + (-0.5)*IT_0144 + 0.5*IT_0220 + (-0.5)
      *IT_0221 + 0.5*IT_0224 + (-0.5)*IT_0225 + 0.5*IT_0228 + (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = IT_0217 + IT_0230;
    const ccomplex_t IT_0232 = IT_0005*IT_0231;
    const ccomplex_t IT_0233 = 1.4142135623731*IT_0232;
    const ccomplex_t IT_0234 = (-0.5)*IT_0233;
    const ccomplex_t IT_0235 = IT_0013 + IT_0017 + IT_0020 + IT_0023 + IT_0024;
    const ccomplex_t IT_0236 = IT_0146*IT_0235;
    const ccomplex_t IT_0237 = IT_0047 + IT_0051 + IT_0054 + IT_0057 + IT_0058;
    const ccomplex_t IT_0238 = IT_0149*IT_0237;
    const ccomplex_t IT_0239 = IT_0072 + IT_0076 + IT_0081 + IT_0082 + IT_0083;
    const ccomplex_t IT_0240 = IT_0157*IT_0239;
    const ccomplex_t IT_0241 = IT_0186 + IT_0188 + IT_0190 + IT_0236 + IT_0238
       + IT_0240;
    const ccomplex_t IT_0242 = IT_0013 + IT_0017 + IT_0019 + IT_0020 + IT_0023
       + IT_0024;
    const ccomplex_t IT_0243 = IT_0035*IT_0242;
    const ccomplex_t IT_0244 = IT_0096*IT_0235;
    const ccomplex_t IT_0245 = IT_0047 + IT_0051 + IT_0053 + IT_0054 + IT_0057
       + IT_0058;
    const ccomplex_t IT_0246 = IT_0042*IT_0245;
    const ccomplex_t IT_0247 = IT_0154*IT_0237;
    const ccomplex_t IT_0248 = IT_0072 + IT_0076 + IT_0080 + IT_0081 + IT_0082
       + IT_0083;
    const ccomplex_t IT_0249 = IT_0067*IT_0248;
    const ccomplex_t IT_0250 = IT_0162*IT_0239;
    const ccomplex_t IT_0251 = -IT_0185 + -IT_0187 + -IT_0189 + 2*IT_0243 + 
      -IT_0244 + 2*IT_0246 + -IT_0247 + 2*IT_0249 + -IT_0250;
    const ccomplex_t IT_0252 = IT_0241 + IT_0251;
    const ccomplex_t IT_0253 = IT_0005*IT_0252;
    const ccomplex_t IT_0254 = 1.4142135623731*IT_0253;
    const ccomplex_t IT_0255 = 0.25*IT_0254;
    return create_ccomplex_return((0 + _Complex_I*0.25)*IT_0168 + (0 +
       _Complex_I*0.25)*IT_0204 + (0 + _Complex_I*0.25)*IT_0234 + (0 +
       _Complex_I*0.25)*IT_0255);
}

