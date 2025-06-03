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
    const ccomplex_t IT_0006 = pow(m_b, 2);
    const ccomplex_t IT_0007 = pow(m_s, 2);
    const ccomplex_t IT_0008 = pow(m_u, 2);
    const ccomplex_t IT_0009 = pow(m_Hp, 2);
    const ccomplex_t IT_0010 = mtylt_C0iC(0, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0008, IT_0008, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0011 = mtylt_C0iC(3, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0008, IT_0008, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0012 = mtylt_C0iC(6, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0008, IT_0008, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0013 = IT_0010 + IT_0011 + IT_0012;
    const ccomplex_t IT_0014 = 0.101321183642338*m_u;
    const ccomplex_t IT_0015 = pow(M_W, -1);
    const ccomplex_t IT_0016 = cos(beta);
    const ccomplex_t IT_0017 = sin(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = sin(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_u*V_ub*e_em
      *IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0026 = 0.666666666666667*IT_0025;
    const ccomplex_t IT_0027 = IT_0022*IT_0024*IT_0026;
    const ccomplex_t IT_0028 = IT_0014*IT_0027;
    const ccomplex_t IT_0029 = IT_0013*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*m_b*V_ub*e_em
      *IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = IT_0026*IT_0031*IT_0033;
    const ccomplex_t IT_0035 = IT_0014*IT_0034;
    const ccomplex_t IT_0036 = IT_0013*IT_0035;
    const ccomplex_t IT_0037 = pow(M_W, 2);
    const ccomplex_t IT_0038 = mtylt_C0iC(3, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0037, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0039 = m_b*IT_0038;
    const ccomplex_t IT_0040 = mtylt_C0iC(6, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0037, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0041 = m_b*IT_0040;
    const ccomplex_t IT_0042 = mtylt_C0iC(12, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0037, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0043 = m_b*IT_0042;
    const ccomplex_t IT_0044 = mtylt_C0iC(15, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0037, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0045 = m_b*IT_0044;
    const ccomplex_t IT_0046 = mtylt_C0iC(18, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0037, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0047 = m_b*IT_0046;
    const ccomplex_t IT_0048 = -IT_0039 + -IT_0041 + -IT_0043 + (-2)*IT_0045 +
       -IT_0047;
    const ccomplex_t IT_0049 = m_s*IT_0040;
    const ccomplex_t IT_0050 = m_s*IT_0044;
    const ccomplex_t IT_0051 = m_s*IT_0046;
    const ccomplex_t IT_0052 = IT_0049 + IT_0050 + IT_0051;
    const ccomplex_t IT_0053 = IT_0048 + IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0015*IT_0020;
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_b*V_ub*e_em
      *IT_0015*IT_0020;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0026*IT_0055*IT_0057;
    const ccomplex_t IT_0059 = 0.101321183642338*IT_0058;
    const ccomplex_t IT_0060 = IT_0053*IT_0059;
    const ccomplex_t IT_0061 = mtylt_C0iC(0, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0008, IT_0037, IT_0008, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0062 = IT_0038 + IT_0040 + IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0015*IT_0020;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0026*IT_0057*IT_0064;
    const ccomplex_t IT_0066 = IT_0014*IT_0065;
    const ccomplex_t IT_0067 = IT_0062*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0015*IT_0020;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0015*IT_0020;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0026*IT_0069*IT_0071;
    const ccomplex_t IT_0073 = 0.101321183642338*IT_0072;
    const ccomplex_t IT_0074 = pow(m_c, 2);
    const ccomplex_t IT_0075 = mtylt_C0iC(6, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0074, IT_0037, IT_0074, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0076 = m_s*IT_0075;
    const ccomplex_t IT_0077 = mtylt_C0iC(15, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0074, IT_0037, IT_0074, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0078 = m_s*IT_0077;
    const ccomplex_t IT_0079 = mtylt_C0iC(18, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0074, IT_0037, IT_0074, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0080 = m_s*IT_0079;
    const ccomplex_t IT_0081 = IT_0076 + IT_0078 + IT_0080;
    const ccomplex_t IT_0082 = mtylt_C0iC(3, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0074, IT_0037, IT_0074, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0083 = m_b*IT_0082;
    const ccomplex_t IT_0084 = m_b*IT_0075;
    const ccomplex_t IT_0085 = m_b*IT_0077;
    const ccomplex_t IT_0086 = m_b*IT_0079;
    const ccomplex_t IT_0087 = mtylt_C0iC(12, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0074, IT_0037, IT_0074, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0088 = m_b*IT_0087;
    const ccomplex_t IT_0089 = -IT_0083 + -IT_0084 + (-2)*IT_0085 + -IT_0086 +
       -IT_0088;
    const ccomplex_t IT_0090 = IT_0081 + IT_0089;
    const ccomplex_t IT_0091 = IT_0073*IT_0090;
    const ccomplex_t IT_0092 = 0.101321183642338*m_c;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = IT_0026*IT_0094*IT_0096;
    const ccomplex_t IT_0098 = IT_0092*IT_0097;
    const ccomplex_t IT_0099 = mtylt_C0iC(0, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0074, IT_0074, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0100 = mtylt_C0iC(3, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0074, IT_0074, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0101 = mtylt_C0iC(6, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0074, IT_0074, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0102 = IT_0099 + IT_0100 + IT_0101;
    const ccomplex_t IT_0103 = IT_0098*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0107 = (-0.5)*IT_0106;
    const ccomplex_t IT_0108 = IT_0026*IT_0105*IT_0107;
    const ccomplex_t IT_0109 = IT_0092*IT_0108;
    const ccomplex_t IT_0110 = IT_0102*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0020;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*V_ub*e_em
      *IT_0020;
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = 0.101321183642338*IT_0115;
    const ccomplex_t IT_0117 = mtylt_C0iC(15, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0008, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0118 = (-2)*IT_0025;
    const ccomplex_t IT_0119 = IT_0117*IT_0118;
    const ccomplex_t IT_0120 = m_s*IT_0119;
    const ccomplex_t IT_0121 = mtylt_C0iC(3, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0008, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0122 = -IT_0025;
    const ccomplex_t IT_0123 = IT_0121*IT_0122;
    const ccomplex_t IT_0124 = m_s*IT_0123;
    const ccomplex_t IT_0125 = mtylt_C0iC(6, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0008, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0126 = IT_0122*IT_0125;
    const ccomplex_t IT_0127 = m_b*IT_0126;
    const ccomplex_t IT_0128 = mtylt_C0iC(18, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0008, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0129 = IT_0118*IT_0128;
    const ccomplex_t IT_0130 = m_b*IT_0129;
    const ccomplex_t IT_0131 = 2*IT_0025;
    const ccomplex_t IT_0132 = IT_0121*IT_0131;
    const ccomplex_t IT_0133 = m_s*IT_0132;
    const ccomplex_t IT_0134 = IT_0125*IT_0131;
    const ccomplex_t IT_0135 = m_s*IT_0134;
    const ccomplex_t IT_0136 = IT_0120 + IT_0124 + IT_0127 + IT_0130 + IT_0133
       + IT_0135;
    const ccomplex_t IT_0137 = IT_0025*IT_0125;
    const ccomplex_t IT_0138 = m_s*IT_0137;
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = IT_0136 + IT_0139;
    const ccomplex_t IT_0141 = IT_0116*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0020;
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0020;
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = IT_0143*IT_0145;
    const ccomplex_t IT_0147 = 0.101321183642338*IT_0146;
    const ccomplex_t IT_0148 = mtylt_C0iC(3, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0074, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0149 = IT_0122*IT_0148;
    const ccomplex_t IT_0150 = m_s*IT_0149;
    const ccomplex_t IT_0151 = mtylt_C0iC(15, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0074, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0152 = IT_0118*IT_0151;
    const ccomplex_t IT_0153 = m_s*IT_0152;
    const ccomplex_t IT_0154 = mtylt_C0iC(6, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0074, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0155 = IT_0122*IT_0154;
    const ccomplex_t IT_0156 = m_b*IT_0155;
    const ccomplex_t IT_0157 = IT_0131*IT_0148;
    const ccomplex_t IT_0158 = m_s*IT_0157;
    const ccomplex_t IT_0159 = IT_0131*IT_0154;
    const ccomplex_t IT_0160 = m_s*IT_0159;
    const ccomplex_t IT_0161 = mtylt_C0iC(18, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0074, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0162 = IT_0118*IT_0161;
    const ccomplex_t IT_0163 = m_b*IT_0162;
    const ccomplex_t IT_0164 = IT_0150 + IT_0153 + IT_0156 + IT_0158 + IT_0160
       + IT_0163;
    const ccomplex_t IT_0165 = IT_0025*IT_0154;
    const ccomplex_t IT_0166 = m_s*IT_0165;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0164 + IT_0167;
    const ccomplex_t IT_0169 = IT_0147*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0020;
    const ccomplex_t IT_0171 = 0.5*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0020;
    const ccomplex_t IT_0173 = 0.5*IT_0172;
    const ccomplex_t IT_0174 = IT_0171*IT_0173;
    const ccomplex_t IT_0175 = 0.101321183642338*IT_0174;
    const ccomplex_t IT_0176 = pow(m_t, 2);
    const ccomplex_t IT_0177 = mtylt_C0iC(6, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0176, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0178 = IT_0131*IT_0177;
    const ccomplex_t IT_0179 = m_s*IT_0178;
    const ccomplex_t IT_0180 = mtylt_C0iC(15, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0176, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0181 = IT_0118*IT_0180;
    const ccomplex_t IT_0182 = m_s*IT_0181;
    const ccomplex_t IT_0183 = mtylt_C0iC(18, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0176, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0184 = IT_0118*IT_0183;
    const ccomplex_t IT_0185 = m_b*IT_0184;
    const ccomplex_t IT_0186 = IT_0122*IT_0177;
    const ccomplex_t IT_0187 = m_b*IT_0186;
    const ccomplex_t IT_0188 = mtylt_C0iC(3, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0176, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0189 = IT_0122*IT_0188;
    const ccomplex_t IT_0190 = m_s*IT_0189;
    const ccomplex_t IT_0191 = IT_0131*IT_0188;
    const ccomplex_t IT_0192 = m_s*IT_0191;
    const ccomplex_t IT_0193 = IT_0179 + IT_0182 + IT_0185 + IT_0187 + IT_0190
       + IT_0192;
    const ccomplex_t IT_0194 = IT_0025*IT_0177;
    const ccomplex_t IT_0195 = m_s*IT_0194;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0193 + IT_0196;
    const ccomplex_t IT_0198 = IT_0175*IT_0197;
    const ccomplex_t IT_0199 = IT_0022*IT_0026*IT_0031;
    const ccomplex_t IT_0200 = 0.101321183642338*IT_0199;
    const ccomplex_t IT_0201 = mtylt_C0iC(15, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0008, IT_0008, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0202 = m_s*IT_0201;
    const ccomplex_t IT_0203 = m_s*IT_0011;
    const ccomplex_t IT_0204 = mtylt_C0iC(12, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0008, IT_0008, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0205 = m_s*IT_0204;
    const ccomplex_t IT_0206 = IT_0202 + IT_0203 + IT_0205;
    const ccomplex_t IT_0207 = m_b*IT_0012;
    const ccomplex_t IT_0208 = m_b*IT_0201;
    const ccomplex_t IT_0209 = mtylt_C0iC(18, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0008, IT_0008, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0210 = m_b*IT_0209;
    const ccomplex_t IT_0211 = m_b*IT_0011;
    const ccomplex_t IT_0212 = m_b*IT_0204;
    const ccomplex_t IT_0213 = -IT_0207 + (-2)*IT_0208 + -IT_0210 + -IT_0211 +
       -IT_0212;
    const ccomplex_t IT_0214 = IT_0206 + IT_0213;
    const ccomplex_t IT_0215 = IT_0200*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1.4142135623731)*m_u*V_ub*e_em
      *IT_0015*IT_0020;
    const ccomplex_t IT_0217 = 0.5*IT_0216;
    const ccomplex_t IT_0218 = IT_0026*IT_0055*IT_0217;
    const ccomplex_t IT_0219 = IT_0014*IT_0218;
    const ccomplex_t IT_0220 = IT_0062*IT_0219;
    const ccomplex_t IT_0221 = IT_0026*IT_0064*IT_0217;
    const ccomplex_t IT_0222 = 0.101321183642338*IT_0221;
    const ccomplex_t IT_0223 = IT_0053*IT_0222;
    const ccomplex_t IT_0224 = IT_0024*IT_0026*IT_0033;
    const ccomplex_t IT_0225 = 0.101321183642338*IT_0224;
    const ccomplex_t IT_0226 = IT_0214*IT_0225;
    const ccomplex_t IT_0227 = IT_0026*IT_0094*IT_0105;
    const ccomplex_t IT_0228 = 0.101321183642338*IT_0227;
    const ccomplex_t IT_0229 = m_s*IT_0100;
    const ccomplex_t IT_0230 = mtylt_C0iC(12, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0074, IT_0074, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0231 = m_s*IT_0230;
    const ccomplex_t IT_0232 = mtylt_C0iC(15, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0074, IT_0074, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0233 = m_s*IT_0232;
    const ccomplex_t IT_0234 = IT_0229 + IT_0231 + IT_0233;
    const ccomplex_t IT_0235 = m_b*IT_0100;
    const ccomplex_t IT_0236 = m_b*IT_0101;
    const ccomplex_t IT_0237 = m_b*IT_0230;
    const ccomplex_t IT_0238 = m_b*IT_0232;
    const ccomplex_t IT_0239 = mtylt_C0iC(18, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0074, IT_0074, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0240 = m_b*IT_0239;
    const ccomplex_t IT_0241 = -IT_0235 + -IT_0236 + -IT_0237 + (-2)*IT_0238 +
       -IT_0240;
    const ccomplex_t IT_0242 = IT_0234 + IT_0241;
    const ccomplex_t IT_0243 = IT_0228*IT_0242;
    const ccomplex_t IT_0244 = mtylt_C0iC(0, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0074, IT_0037, IT_0074, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0245 = IT_0075 + IT_0082 + IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0015*IT_0020;
    const ccomplex_t IT_0247 = 0.5*IT_0246;
    const ccomplex_t IT_0248 = IT_0026*IT_0071*IT_0247;
    const ccomplex_t IT_0249 = IT_0092*IT_0248;
    const ccomplex_t IT_0250 = IT_0245*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0015*IT_0020;
    const ccomplex_t IT_0252 = 0.5*IT_0251;
    const ccomplex_t IT_0253 = IT_0026*IT_0247*IT_0252;
    const ccomplex_t IT_0254 = 0.101321183642338*IT_0253;
    const ccomplex_t IT_0255 = IT_0090*IT_0254;
    const ccomplex_t IT_0256 = IT_0026*IT_0096*IT_0107;
    const ccomplex_t IT_0257 = 0.101321183642338*IT_0256;
    const ccomplex_t IT_0258 = IT_0242*IT_0257;
    const ccomplex_t IT_0259 = IT_0026*IT_0069*IT_0252;
    const ccomplex_t IT_0260 = IT_0092*IT_0259;
    const ccomplex_t IT_0261 = IT_0245*IT_0260;
    const ccomplex_t IT_0262 = 0.101321183642338*m_t;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0264 = 0.5*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0266 = (-0.5)*IT_0265;
    const ccomplex_t IT_0267 = IT_0026*IT_0264*IT_0266;
    const ccomplex_t IT_0268 = IT_0262*IT_0267;
    const ccomplex_t IT_0269 = mtylt_C0iC(0, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0176, IT_0176, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0270 = mtylt_C0iC(3, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0176, IT_0176, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0271 = mtylt_C0iC(6, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0176, IT_0176, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0272 = IT_0269 + IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = IT_0268*IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0275 = 0.5*IT_0274;
    const ccomplex_t IT_0276 = IT_0026*IT_0264*IT_0275;
    const ccomplex_t IT_0277 = 0.101321183642338*IT_0276;
    const ccomplex_t IT_0278 = m_s*IT_0270;
    const ccomplex_t IT_0279 = mtylt_C0iC(12, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0176, IT_0176, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0280 = m_s*IT_0279;
    const ccomplex_t IT_0281 = mtylt_C0iC(15, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0176, IT_0176, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0282 = m_s*IT_0281;
    const ccomplex_t IT_0283 = IT_0278 + IT_0280 + IT_0282;
    const ccomplex_t IT_0284 = m_b*IT_0270;
    const ccomplex_t IT_0285 = m_b*IT_0279;
    const ccomplex_t IT_0286 = m_b*IT_0281;
    const ccomplex_t IT_0287 = m_b*IT_0271;
    const ccomplex_t IT_0288 = mtylt_C0iC(18, IT_0006 + IT_0007 + (-2)*s_12,
       IT_0007, IT_0006, IT_0176, IT_0176, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0289 = m_b*IT_0288;
    const ccomplex_t IT_0290 = -IT_0284 + -IT_0285 + (-2)*IT_0286 + -IT_0287 +
       -IT_0289;
    const ccomplex_t IT_0291 = IT_0283 + IT_0290;
    const ccomplex_t IT_0292 = IT_0277*IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0015*IT_0020;
    const ccomplex_t IT_0294 = 0.5*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0015*IT_0020;
    const ccomplex_t IT_0296 = (-0.5)*IT_0295;
    const ccomplex_t IT_0297 = IT_0026*IT_0294*IT_0296;
    const ccomplex_t IT_0298 = IT_0262*IT_0297;
    const ccomplex_t IT_0299 = mtylt_C0iC(3, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0176, IT_0037, IT_0176, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0300 = mtylt_C0iC(6, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0176, IT_0037, IT_0176, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0301 = mtylt_C0iC(0, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0176, IT_0037, IT_0176, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0302 = IT_0299 + IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = IT_0298*IT_0302;
    const ccomplex_t IT_0304 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0015*IT_0020;
    const ccomplex_t IT_0305 = 0.5*IT_0304;
    const ccomplex_t IT_0306 = IT_0026*IT_0294*IT_0305;
    const ccomplex_t IT_0307 = 0.101321183642338*IT_0306;
    const ccomplex_t IT_0308 = m_b*IT_0299;
    const ccomplex_t IT_0309 = m_b*IT_0300;
    const ccomplex_t IT_0310 = mtylt_C0iC(12, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0176, IT_0037, IT_0176, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0311 = m_b*IT_0310;
    const ccomplex_t IT_0312 = mtylt_C0iC(15, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0176, IT_0037, IT_0176, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0313 = m_b*IT_0312;
    const ccomplex_t IT_0314 = mtylt_C0iC(18, IT_0006, IT_0007, IT_0006 +
       IT_0007 + (-2)*s_12, IT_0176, IT_0037, IT_0176, getIntegrationParameters(
      )->reg_int);
    const ccomplex_t IT_0315 = m_b*IT_0314;
    const ccomplex_t IT_0316 = -IT_0308 + -IT_0309 + -IT_0311 + (-2)*IT_0313 +
       -IT_0315;
    const ccomplex_t IT_0317 = m_s*IT_0300;
    const ccomplex_t IT_0318 = m_s*IT_0312;
    const ccomplex_t IT_0319 = m_s*IT_0314;
    const ccomplex_t IT_0320 = IT_0317 + IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = IT_0316 + IT_0320;
    const ccomplex_t IT_0322 = IT_0307*IT_0321;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0015*IT_0016*IT_0018*IT_0020;
    const ccomplex_t IT_0324 = (-0.5)*IT_0323;
    const ccomplex_t IT_0325 = IT_0026*IT_0266*IT_0324;
    const ccomplex_t IT_0326 = 0.101321183642338*IT_0325;
    const ccomplex_t IT_0327 = IT_0291*IT_0326;
    const ccomplex_t IT_0328 = IT_0026*IT_0275*IT_0324;
    const ccomplex_t IT_0329 = IT_0262*IT_0328;
    const ccomplex_t IT_0330 = IT_0272*IT_0329;
    const ccomplex_t IT_0331 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0015*IT_0020;
    const ccomplex_t IT_0332 = (-0.5)*IT_0331;
    const ccomplex_t IT_0333 = IT_0026*IT_0296*IT_0332;
    const ccomplex_t IT_0334 = 0.101321183642338*IT_0333;
    const ccomplex_t IT_0335 = IT_0321*IT_0334;
    const ccomplex_t IT_0336 = IT_0026*IT_0305*IT_0332;
    const ccomplex_t IT_0337 = IT_0262*IT_0336;
    const ccomplex_t IT_0338 = IT_0302*IT_0337;
    const ccomplex_t IT_0339 = IT_0029 + IT_0036 + IT_0060 + IT_0067 + IT_0091
       + IT_0103 + IT_0110 + IT_0141 + IT_0169 + IT_0198 + IT_0215 + IT_0220 +
       IT_0223 + IT_0226 + IT_0243 + IT_0250 + IT_0255 + IT_0258 + IT_0261 +
       IT_0273 + IT_0292 + IT_0303 + IT_0322 + IT_0327 + IT_0330 + IT_0335 +
       IT_0338;
    const ccomplex_t IT_0340 = IT_0026*IT_0112*IT_0114;
    const ccomplex_t IT_0341 = 0.101321183642338*IT_0340;
    const ccomplex_t IT_0342 = m_s*IT_0038;
    const ccomplex_t IT_0343 = IT_0049 + IT_0050 + IT_0051 + IT_0342;
    const ccomplex_t IT_0344 = IT_0048 + IT_0343;
    const ccomplex_t IT_0345 = IT_0341*IT_0344;
    const ccomplex_t IT_0346 = (0 + _Complex_I*1)*M_W*e_em;
    const ccomplex_t IT_0347 = IT_0055*IT_0114*IT_0346;
    const ccomplex_t IT_0348 = 0.101321183642338*IT_0347;
    const ccomplex_t IT_0349 = IT_0125*IT_0348;
    const ccomplex_t IT_0350 = IT_0071*IT_0145*IT_0346;
    const ccomplex_t IT_0351 = 0.101321183642338*IT_0350;
    const ccomplex_t IT_0352 = IT_0154*IT_0351;
    const ccomplex_t IT_0353 = IT_0173*IT_0296*IT_0346;
    const ccomplex_t IT_0354 = 0.101321183642338*IT_0353;
    const ccomplex_t IT_0355 = IT_0177*IT_0354;
    const ccomplex_t IT_0356 = IT_0022*IT_0024;
    const ccomplex_t IT_0357 = IT_0014*IT_0356;
    const ccomplex_t IT_0358 = mtylt_C0iC(0, IT_0006, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0007, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0359 = IT_0122*IT_0358;
    const ccomplex_t IT_0360 = IT_0118*IT_0358;
    const ccomplex_t IT_0361 = mtylt_C0iC(3, IT_0006, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0007, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0362 = IT_0118*IT_0361;
    const ccomplex_t IT_0363 = -IT_0360 + -IT_0362;
    const ccomplex_t IT_0364 = IT_0359 + IT_0363;
    const ccomplex_t IT_0365 = IT_0357*IT_0364;
    const ccomplex_t IT_0366 = IT_0022*IT_0031;
    const ccomplex_t IT_0367 = 0.101321183642338*IT_0366;
    const ccomplex_t IT_0368 = mtylt_C0iC(6, IT_0006, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0007, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0369 = IT_0122*IT_0368;
    const ccomplex_t IT_0370 = m_s*IT_0369;
    const ccomplex_t IT_0371 = IT_0122*IT_0361;
    const ccomplex_t IT_0372 = m_b*IT_0371;
    const ccomplex_t IT_0373 = IT_0370 + IT_0372;
    const ccomplex_t IT_0374 = m_b*IT_0362;
    const ccomplex_t IT_0375 = mtylt_C0iC(12, IT_0006, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0007, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0376 = IT_0118*IT_0375;
    const ccomplex_t IT_0377 = m_b*IT_0376;
    const ccomplex_t IT_0378 = IT_0118*IT_0368;
    const ccomplex_t IT_0379 = m_s*IT_0378;
    const ccomplex_t IT_0380 = mtylt_C0iC(15, IT_0006, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0007, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0381 = IT_0118*IT_0380;
    const ccomplex_t IT_0382 = m_s*IT_0381;
    const ccomplex_t IT_0383 = -IT_0374 + -IT_0377 + -IT_0379 + -IT_0382;
    const ccomplex_t IT_0384 = IT_0373 + IT_0383;
    const ccomplex_t IT_0385 = IT_0367*IT_0384;
    const ccomplex_t IT_0386 = IT_0094*IT_0096;
    const ccomplex_t IT_0387 = IT_0092*IT_0386;
    const ccomplex_t IT_0388 = mtylt_C0iC(0, IT_0006, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0007, IT_0074, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0389 = IT_0122*IT_0388;
    const ccomplex_t IT_0390 = IT_0118*IT_0388;
    const ccomplex_t IT_0391 = mtylt_C0iC(3, IT_0006, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0007, IT_0074, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0392 = IT_0118*IT_0391;
    const ccomplex_t IT_0393 = -IT_0390 + -IT_0392;
    const ccomplex_t IT_0394 = IT_0389 + IT_0393;
    const ccomplex_t IT_0395 = IT_0387*IT_0394;
    const ccomplex_t IT_0396 = IT_0094*IT_0105;
    const ccomplex_t IT_0397 = 0.101321183642338*IT_0396;
    const ccomplex_t IT_0398 = mtylt_C0iC(6, IT_0006, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0007, IT_0074, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0399 = IT_0122*IT_0398;
    const ccomplex_t IT_0400 = m_s*IT_0399;
    const ccomplex_t IT_0401 = IT_0122*IT_0391;
    const ccomplex_t IT_0402 = m_b*IT_0401;
    const ccomplex_t IT_0403 = IT_0400 + IT_0402;
    const ccomplex_t IT_0404 = m_b*IT_0392;
    const ccomplex_t IT_0405 = mtylt_C0iC(12, IT_0006, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0007, IT_0074, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0406 = IT_0118*IT_0405;
    const ccomplex_t IT_0407 = m_b*IT_0406;
    const ccomplex_t IT_0408 = IT_0118*IT_0398;
    const ccomplex_t IT_0409 = m_s*IT_0408;
    const ccomplex_t IT_0410 = mtylt_C0iC(15, IT_0006, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0007, IT_0074, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0411 = IT_0118*IT_0410;
    const ccomplex_t IT_0412 = m_s*IT_0411;
    const ccomplex_t IT_0413 = -IT_0404 + -IT_0407 + -IT_0409 + -IT_0412;
    const ccomplex_t IT_0414 = IT_0403 + IT_0413;
    const ccomplex_t IT_0415 = IT_0397*IT_0414;
    const ccomplex_t IT_0416 = IT_0264*IT_0266;
    const ccomplex_t IT_0417 = IT_0262*IT_0416;
    const ccomplex_t IT_0418 = mtylt_C0iC(0, IT_0006, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0007, IT_0176, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0419 = IT_0118*IT_0418;
    const ccomplex_t IT_0420 = mtylt_C0iC(3, IT_0006, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0007, IT_0176, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0421 = IT_0118*IT_0420;
    const ccomplex_t IT_0422 = -IT_0419 + -IT_0421;
    const ccomplex_t IT_0423 = IT_0122*IT_0418;
    const ccomplex_t IT_0424 = IT_0422 + IT_0423;
    const ccomplex_t IT_0425 = IT_0417*IT_0424;
    const ccomplex_t IT_0426 = IT_0264*IT_0275;
    const ccomplex_t IT_0427 = 0.101321183642338*IT_0426;
    const ccomplex_t IT_0428 = mtylt_C0iC(6, IT_0006, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0007, IT_0176, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0429 = IT_0122*IT_0428;
    const ccomplex_t IT_0430 = m_s*IT_0429;
    const ccomplex_t IT_0431 = IT_0122*IT_0420;
    const ccomplex_t IT_0432 = m_b*IT_0431;
    const ccomplex_t IT_0433 = IT_0430 + IT_0432;
    const ccomplex_t IT_0434 = IT_0118*IT_0428;
    const ccomplex_t IT_0435 = m_s*IT_0434;
    const ccomplex_t IT_0436 = m_b*IT_0421;
    const ccomplex_t IT_0437 = mtylt_C0iC(12, IT_0006, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0007, IT_0176, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0438 = IT_0118*IT_0437;
    const ccomplex_t IT_0439 = m_b*IT_0438;
    const ccomplex_t IT_0440 = mtylt_C0iC(15, IT_0006, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0007, IT_0176, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0441 = IT_0118*IT_0440;
    const ccomplex_t IT_0442 = m_s*IT_0441;
    const ccomplex_t IT_0443 = -IT_0435 + -IT_0436 + -IT_0439 + -IT_0442;
    const ccomplex_t IT_0444 = IT_0433 + IT_0443;
    const ccomplex_t IT_0445 = IT_0427*IT_0444;
    const ccomplex_t IT_0446 = IT_0024*IT_0033;
    const ccomplex_t IT_0447 = 0.101321183642338*IT_0446;
    const ccomplex_t IT_0448 = IT_0384*IT_0447;
    const ccomplex_t IT_0449 = IT_0031*IT_0033;
    const ccomplex_t IT_0450 = IT_0014*IT_0449;
    const ccomplex_t IT_0451 = IT_0364*IT_0450;
    const ccomplex_t IT_0452 = IT_0096*IT_0107;
    const ccomplex_t IT_0453 = 0.101321183642338*IT_0452;
    const ccomplex_t IT_0454 = IT_0414*IT_0453;
    const ccomplex_t IT_0455 = IT_0105*IT_0107;
    const ccomplex_t IT_0456 = IT_0092*IT_0455;
    const ccomplex_t IT_0457 = IT_0394*IT_0456;
    const ccomplex_t IT_0458 = IT_0266*IT_0324;
    const ccomplex_t IT_0459 = 0.101321183642338*IT_0458;
    const ccomplex_t IT_0460 = IT_0444*IT_0459;
    const ccomplex_t IT_0461 = IT_0275*IT_0324;
    const ccomplex_t IT_0462 = IT_0262*IT_0461;
    const ccomplex_t IT_0463 = IT_0424*IT_0462;
    const ccomplex_t IT_0464 = IT_0026*IT_0143*IT_0145;
    const ccomplex_t IT_0465 = 0.101321183642338*IT_0464;
    const ccomplex_t IT_0466 = m_s*IT_0082;
    const ccomplex_t IT_0467 = IT_0076 + IT_0078 + IT_0080 + IT_0466;
    const ccomplex_t IT_0468 = IT_0089 + IT_0467;
    const ccomplex_t IT_0469 = IT_0465*IT_0468;
    const ccomplex_t IT_0470 = IT_0026*IT_0171*IT_0173;
    const ccomplex_t IT_0471 = 0.101321183642338*IT_0470;
    const ccomplex_t IT_0472 = m_s*IT_0299;
    const ccomplex_t IT_0473 = IT_0317 + IT_0318 + IT_0319 + IT_0472;
    const ccomplex_t IT_0474 = IT_0316 + IT_0473;
    const ccomplex_t IT_0475 = IT_0471*IT_0474;
    const ccomplex_t IT_0476 = IT_0055*IT_0217;
    const ccomplex_t IT_0477 = IT_0014*IT_0476;
    const ccomplex_t IT_0478 = mtylt_C0iC(0, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0008, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0479 = IT_0118*IT_0478;
    const ccomplex_t IT_0480 = IT_0118*IT_0125;
    const ccomplex_t IT_0481 = -IT_0479 + -IT_0480;
    const ccomplex_t IT_0482 = IT_0122*IT_0478;
    const ccomplex_t IT_0483 = IT_0481 + IT_0482;
    const ccomplex_t IT_0484 = IT_0477*IT_0483;
    const ccomplex_t IT_0485 = IT_0064*IT_0217;
    const ccomplex_t IT_0486 = 0.101321183642338*IT_0485;
    const ccomplex_t IT_0487 = IT_0124 + IT_0127;
    const ccomplex_t IT_0488 = IT_0118*IT_0121;
    const ccomplex_t IT_0489 = m_s*IT_0488;
    const ccomplex_t IT_0490 = m_b*IT_0480;
    const ccomplex_t IT_0491 = -IT_0120 + -IT_0130 + -IT_0489 + -IT_0490;
    const ccomplex_t IT_0492 = IT_0487 + IT_0491;
    const ccomplex_t IT_0493 = IT_0486*IT_0492;
    const ccomplex_t IT_0494 = mtylt_C0iC(0, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0074, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0495 = IT_0122*IT_0494;
    const ccomplex_t IT_0496 = IT_0118*IT_0494;
    const ccomplex_t IT_0497 = IT_0118*IT_0154;
    const ccomplex_t IT_0498 = -IT_0496 + -IT_0497;
    const ccomplex_t IT_0499 = IT_0495 + IT_0498;
    const ccomplex_t IT_0500 = IT_0071*IT_0247;
    const ccomplex_t IT_0501 = IT_0092*IT_0500;
    const ccomplex_t IT_0502 = IT_0499*IT_0501;
    const ccomplex_t IT_0503 = IT_0247*IT_0252;
    const ccomplex_t IT_0504 = 0.101321183642338*IT_0503;
    const ccomplex_t IT_0505 = IT_0150 + IT_0156;
    const ccomplex_t IT_0506 = m_b*IT_0497;
    const ccomplex_t IT_0507 = IT_0118*IT_0148;
    const ccomplex_t IT_0508 = m_s*IT_0507;
    const ccomplex_t IT_0509 = -IT_0153 + -IT_0163 + -IT_0506 + -IT_0508;
    const ccomplex_t IT_0510 = IT_0505 + IT_0509;
    const ccomplex_t IT_0511 = IT_0504*IT_0510;
    const ccomplex_t IT_0512 = IT_0294*IT_0296;
    const ccomplex_t IT_0513 = IT_0262*IT_0512;
    const ccomplex_t IT_0514 = mtylt_C0iC(0, IT_0007, IT_0006 + IT_0007 + (-2)
      *s_12, IT_0006, IT_0176, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0515 = IT_0122*IT_0514;
    const ccomplex_t IT_0516 = IT_0118*IT_0514;
    const ccomplex_t IT_0517 = IT_0118*IT_0177;
    const ccomplex_t IT_0518 = -IT_0516 + -IT_0517;
    const ccomplex_t IT_0519 = IT_0515 + IT_0518;
    const ccomplex_t IT_0520 = IT_0513*IT_0519;
    const ccomplex_t IT_0521 = IT_0294*IT_0305;
    const ccomplex_t IT_0522 = 0.101321183642338*IT_0521;
    const ccomplex_t IT_0523 = IT_0187 + IT_0190;
    const ccomplex_t IT_0524 = IT_0118*IT_0188;
    const ccomplex_t IT_0525 = m_s*IT_0524;
    const ccomplex_t IT_0526 = m_b*IT_0517;
    const ccomplex_t IT_0527 = -IT_0182 + -IT_0185 + -IT_0525 + -IT_0526;
    const ccomplex_t IT_0528 = IT_0523 + IT_0527;
    const ccomplex_t IT_0529 = IT_0522*IT_0528;
    const ccomplex_t IT_0530 = IT_0055*IT_0057;
    const ccomplex_t IT_0531 = 0.101321183642338*IT_0530;
    const ccomplex_t IT_0532 = IT_0492*IT_0531;
    const ccomplex_t IT_0533 = IT_0057*IT_0064;
    const ccomplex_t IT_0534 = IT_0014*IT_0533;
    const ccomplex_t IT_0535 = IT_0483*IT_0534;
    const ccomplex_t IT_0536 = IT_0069*IT_0071;
    const ccomplex_t IT_0537 = 0.101321183642338*IT_0536;
    const ccomplex_t IT_0538 = IT_0510*IT_0537;
    const ccomplex_t IT_0539 = IT_0069*IT_0252;
    const ccomplex_t IT_0540 = IT_0092*IT_0539;
    const ccomplex_t IT_0541 = IT_0499*IT_0540;
    const ccomplex_t IT_0542 = IT_0296*IT_0332;
    const ccomplex_t IT_0543 = 0.101321183642338*IT_0542;
    const ccomplex_t IT_0544 = IT_0528*IT_0543;
    const ccomplex_t IT_0545 = IT_0305*IT_0332;
    const ccomplex_t IT_0546 = IT_0262*IT_0545;
    const ccomplex_t IT_0547 = IT_0519*IT_0546;
    const ccomplex_t IT_0548 = 2*IT_0345 + -IT_0349 + -IT_0352 + -IT_0355 +
       0.5*IT_0365 + (-0.5)*IT_0385 + 0.5*IT_0395 + (-0.5)*IT_0415 + 0.5*IT_0425
       + (-0.5)*IT_0445 + (-0.5)*IT_0448 + 0.5*IT_0451 + (-0.5)*IT_0454 + 0.5
      *IT_0457 + (-0.5)*IT_0460 + 0.5*IT_0463 + 2*IT_0469 + 2*IT_0475 + 0.5
      *IT_0484 + (-0.5)*IT_0493 + 0.5*IT_0502 + (-0.5)*IT_0511 + 0.5*IT_0520 + (
      -0.5)*IT_0529 + (-0.5)*IT_0532 + 0.5*IT_0535 + (-0.5)*IT_0538 + 0.5
      *IT_0541 + (-0.5)*IT_0544 + 0.5*IT_0547;
    const ccomplex_t IT_0549 = IT_0339 + IT_0548;
    const ccomplex_t IT_0550 = IT_0005*IT_0549;
    const ccomplex_t IT_0551 = 1.4142135623731*IT_0550;
    const ccomplex_t IT_0552 = 0.25*IT_0551;
    const ccomplex_t IT_0553 = IT_0011*IT_0028;
    const ccomplex_t IT_0554 = IT_0011*IT_0035;
    const ccomplex_t IT_0555 = IT_0049 + IT_0051;
    const ccomplex_t IT_0556 = -IT_0041 + -IT_0045 + -IT_0047;
    const ccomplex_t IT_0557 = IT_0555 + IT_0556;
    const ccomplex_t IT_0558 = IT_0059*IT_0557;
    const ccomplex_t IT_0559 = IT_0100*IT_0109;
    const ccomplex_t IT_0560 = IT_0268*IT_0270;
    const ccomplex_t IT_0561 = IT_0057*IT_0112*IT_0346;
    const ccomplex_t IT_0562 = 0.101321183642338*IT_0561;
    const ccomplex_t IT_0563 = IT_0121*IT_0562;
    const ccomplex_t IT_0564 = IT_0069*IT_0143*IT_0346;
    const ccomplex_t IT_0565 = 0.101321183642338*IT_0564;
    const ccomplex_t IT_0566 = IT_0148*IT_0565;
    const ccomplex_t IT_0567 = IT_0171*IT_0332*IT_0346;
    const ccomplex_t IT_0568 = 0.101321183642338*IT_0567;
    const ccomplex_t IT_0569 = IT_0188*IT_0568;
    const ccomplex_t IT_0570 = IT_0040*IT_0066;
    const ccomplex_t IT_0571 = IT_0098*IT_0100;
    const ccomplex_t IT_0572 = IT_0075*IT_0249;
    const ccomplex_t IT_0573 = m_b*IT_0123;
    const ccomplex_t IT_0574 = IT_0124 + IT_0573;
    const ccomplex_t IT_0575 = m_b*IT_0137;
    const ccomplex_t IT_0576 = m_b*IT_0119;
    const ccomplex_t IT_0577 = mtylt_C0iC(12, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0008, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0578 = IT_0118*IT_0577;
    const ccomplex_t IT_0579 = m_s*IT_0578;
    const ccomplex_t IT_0580 = 2*IT_0127 + -IT_0489 + -IT_0490 + -IT_0575 + 
      -IT_0576 + -IT_0579;
    const ccomplex_t IT_0581 = IT_0574 + IT_0580;
    const ccomplex_t IT_0582 = IT_0116*IT_0581;
    const ccomplex_t IT_0583 = m_b*IT_0149;
    const ccomplex_t IT_0584 = IT_0150 + IT_0583;
    const ccomplex_t IT_0585 = mtylt_C0iC(12, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0074, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0586 = IT_0118*IT_0585;
    const ccomplex_t IT_0587 = m_s*IT_0586;
    const ccomplex_t IT_0588 = m_b*IT_0165;
    const ccomplex_t IT_0589 = m_b*IT_0152;
    const ccomplex_t IT_0590 = 2*IT_0156 + -IT_0506 + -IT_0508 + -IT_0587 + 
      -IT_0588 + -IT_0589;
    const ccomplex_t IT_0591 = IT_0584 + IT_0590;
    const ccomplex_t IT_0592 = IT_0147*IT_0591;
    const ccomplex_t IT_0593 = m_b*IT_0189;
    const ccomplex_t IT_0594 = IT_0190 + IT_0593;
    const ccomplex_t IT_0595 = m_b*IT_0181;
    const ccomplex_t IT_0596 = mtylt_C0iC(12, IT_0007, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0006, IT_0176, IT_0037, IT_0037, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0597 = IT_0118*IT_0596;
    const ccomplex_t IT_0598 = m_s*IT_0597;
    const ccomplex_t IT_0599 = m_b*IT_0194;
    const ccomplex_t IT_0600 = 2*IT_0187 + -IT_0525 + -IT_0526 + -IT_0595 + 
      -IT_0598 + -IT_0599;
    const ccomplex_t IT_0601 = IT_0594 + IT_0600;
    const ccomplex_t IT_0602 = IT_0175*IT_0601;
    const ccomplex_t IT_0603 = IT_0203 + IT_0205;
    const ccomplex_t IT_0604 = -IT_0208 + -IT_0211 + -IT_0212;
    const ccomplex_t IT_0605 = IT_0603 + IT_0604;
    const ccomplex_t IT_0606 = IT_0200*IT_0605;
    const ccomplex_t IT_0607 = IT_0040*IT_0219;
    const ccomplex_t IT_0608 = IT_0222*IT_0557;
    const ccomplex_t IT_0609 = IT_0225*IT_0605;
    const ccomplex_t IT_0610 = IT_0229 + IT_0231;
    const ccomplex_t IT_0611 = -IT_0235 + -IT_0237 + -IT_0238;
    const ccomplex_t IT_0612 = IT_0610 + IT_0611;
    const ccomplex_t IT_0613 = IT_0228*IT_0612;
    const ccomplex_t IT_0614 = IT_0076 + IT_0080;
    const ccomplex_t IT_0615 = -IT_0084 + -IT_0085 + -IT_0086;
    const ccomplex_t IT_0616 = IT_0614 + IT_0615;
    const ccomplex_t IT_0617 = IT_0254*IT_0616;
    const ccomplex_t IT_0618 = IT_0257*IT_0612;
    const ccomplex_t IT_0619 = IT_0073*IT_0616;
    const ccomplex_t IT_0620 = IT_0075*IT_0260;
    const ccomplex_t IT_0621 = IT_0278 + IT_0280;
    const ccomplex_t IT_0622 = -IT_0284 + -IT_0285 + -IT_0286;
    const ccomplex_t IT_0623 = IT_0621 + IT_0622;
    const ccomplex_t IT_0624 = IT_0277*IT_0623;
    const ccomplex_t IT_0625 = IT_0298*IT_0300;
    const ccomplex_t IT_0626 = IT_0317 + IT_0319;
    const ccomplex_t IT_0627 = -IT_0309 + -IT_0313 + -IT_0315;
    const ccomplex_t IT_0628 = IT_0626 + IT_0627;
    const ccomplex_t IT_0629 = IT_0307*IT_0628;
    const ccomplex_t IT_0630 = IT_0326*IT_0623;
    const ccomplex_t IT_0631 = IT_0270*IT_0329;
    const ccomplex_t IT_0632 = IT_0334*IT_0628;
    const ccomplex_t IT_0633 = IT_0300*IT_0337;
    const ccomplex_t IT_0634 = IT_0553 + IT_0554 + IT_0558 + IT_0559 + IT_0560
       + IT_0563 + IT_0566 + IT_0569 + IT_0570 + IT_0571 + IT_0572 + IT_0582 +
       IT_0592 + IT_0602 + IT_0606 + IT_0607 + IT_0608 + IT_0609 + IT_0613 +
       IT_0617 + IT_0618 + IT_0619 + IT_0620 + IT_0624 + IT_0625 + IT_0629 +
       IT_0630 + IT_0631 + IT_0632 + IT_0633;
    const ccomplex_t IT_0635 = IT_0423 + IT_0434;
    const ccomplex_t IT_0636 = IT_0417*IT_0635;
    const ccomplex_t IT_0637 = IT_0462*IT_0635;
    const ccomplex_t IT_0638 = IT_0482 + IT_0488;
    const ccomplex_t IT_0639 = IT_0477*IT_0638;
    const ccomplex_t IT_0640 = IT_0359 + IT_0378;
    const ccomplex_t IT_0641 = IT_0357*IT_0640;
    const ccomplex_t IT_0642 = mtylt_C0iC(18, IT_0006, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0007, IT_0008, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0643 = IT_0118*IT_0642;
    const ccomplex_t IT_0644 = m_s*IT_0643;
    const ccomplex_t IT_0645 = m_b*IT_0381;
    const ccomplex_t IT_0646 = IT_0370 + IT_0372 + IT_0644 + IT_0645;
    const ccomplex_t IT_0647 = IT_0367*IT_0646;
    const ccomplex_t IT_0648 = IT_0389 + IT_0408;
    const ccomplex_t IT_0649 = IT_0387*IT_0648;
    const ccomplex_t IT_0650 = mtylt_C0iC(18, IT_0006, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0007, IT_0074, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0651 = IT_0118*IT_0650;
    const ccomplex_t IT_0652 = m_s*IT_0651;
    const ccomplex_t IT_0653 = m_b*IT_0411;
    const ccomplex_t IT_0654 = IT_0400 + IT_0402 + IT_0652 + IT_0653;
    const ccomplex_t IT_0655 = IT_0397*IT_0654;
    const ccomplex_t IT_0656 = mtylt_C0iC(18, IT_0006, IT_0006 + IT_0007 + (-2
      )*s_12, IT_0007, IT_0176, IT_0009, IT_0009, getIntegrationParameters()
      ->reg_int);
    const ccomplex_t IT_0657 = IT_0118*IT_0656;
    const ccomplex_t IT_0658 = m_s*IT_0657;
    const ccomplex_t IT_0659 = m_b*IT_0441;
    const ccomplex_t IT_0660 = IT_0430 + IT_0432 + IT_0658 + IT_0659;
    const ccomplex_t IT_0661 = IT_0427*IT_0660;
    const ccomplex_t IT_0662 = IT_0447*IT_0646;
    const ccomplex_t IT_0663 = IT_0450*IT_0640;
    const ccomplex_t IT_0664 = IT_0453*IT_0654;
    const ccomplex_t IT_0665 = IT_0456*IT_0648;
    const ccomplex_t IT_0666 = IT_0459*IT_0660;
    const ccomplex_t IT_0667 = -IT_0039 + -IT_0041 + -IT_0045 + -IT_0047;
    const ccomplex_t IT_0668 = IT_0555 + IT_0667;
    const ccomplex_t IT_0669 = IT_0341*IT_0668;
    const ccomplex_t IT_0670 = -IT_0083 + -IT_0084 + -IT_0085 + -IT_0086;
    const ccomplex_t IT_0671 = IT_0614 + IT_0670;
    const ccomplex_t IT_0672 = IT_0465*IT_0671;
    const ccomplex_t IT_0673 = -IT_0308 + -IT_0309 + -IT_0313 + -IT_0315;
    const ccomplex_t IT_0674 = IT_0626 + IT_0673;
    const ccomplex_t IT_0675 = IT_0471*IT_0674;
    const ccomplex_t IT_0676 = IT_0124 + IT_0127 + IT_0576 + IT_0579;
    const ccomplex_t IT_0677 = IT_0486*IT_0676;
    const ccomplex_t IT_0678 = IT_0495 + IT_0507;
    const ccomplex_t IT_0679 = IT_0501*IT_0678;
    const ccomplex_t IT_0680 = IT_0150 + IT_0156 + IT_0587 + IT_0589;
    const ccomplex_t IT_0681 = IT_0504*IT_0680;
    const ccomplex_t IT_0682 = IT_0515 + IT_0524;
    const ccomplex_t IT_0683 = IT_0513*IT_0682;
    const ccomplex_t IT_0684 = IT_0187 + IT_0190 + IT_0595 + IT_0598;
    const ccomplex_t IT_0685 = IT_0522*IT_0684;
    const ccomplex_t IT_0686 = IT_0531*IT_0676;
    const ccomplex_t IT_0687 = IT_0534*IT_0638;
    const ccomplex_t IT_0688 = IT_0537*IT_0680;
    const ccomplex_t IT_0689 = IT_0540*IT_0678;
    const ccomplex_t IT_0690 = IT_0543*IT_0684;
    const ccomplex_t IT_0691 = IT_0546*IT_0682;
    const ccomplex_t IT_0692 = 0.5*IT_0636 + 0.5*IT_0637 + 0.5*IT_0639 + 0.5
      *IT_0641 + (-0.5)*IT_0647 + 0.5*IT_0649 + (-0.5)*IT_0655 + (-0.5)*IT_0661 
      + (-0.5)*IT_0662 + 0.5*IT_0663 + (-0.5)*IT_0664 + 0.5*IT_0665 + (-0.5)
      *IT_0666 + 2*IT_0669 + 2*IT_0672 + 2*IT_0675 + (-0.5)*IT_0677 + 0.5
      *IT_0679 + (-0.5)*IT_0681 + 0.5*IT_0683 + (-0.5)*IT_0685 + (-0.5)*IT_0686 
      + 0.5*IT_0687 + (-0.5)*IT_0688 + 0.5*IT_0689 + (-0.5)*IT_0690 + 0.5*IT_0691;
    const ccomplex_t IT_0693 = IT_0634 + IT_0692;
    const ccomplex_t IT_0694 = IT_0005*IT_0693;
    const ccomplex_t IT_0695 = 1.4142135623731*IT_0694;
    const ccomplex_t IT_0696 = (-0.25)*IT_0695;
    const ccomplex_t IT_0697 = IT_0229 + IT_0231 + IT_0233 + IT_0235 + IT_0236
       + IT_0237 + IT_0240;
    const ccomplex_t IT_0698 = 2*IT_0238;
    const ccomplex_t IT_0699 = IT_0697 + IT_0698;
    const ccomplex_t IT_0700 = IT_0228*IT_0699;
    const ccomplex_t IT_0701 = IT_0120 + IT_0124 + IT_0133 + IT_0135;
    const ccomplex_t IT_0702 = -IT_0127 + -IT_0130 + -IT_0138;
    const ccomplex_t IT_0703 = IT_0701 + IT_0702;
    const ccomplex_t IT_0704 = IT_0116*IT_0703;
    const ccomplex_t IT_0705 = IT_0150 + IT_0153 + IT_0158 + IT_0160;
    const ccomplex_t IT_0706 = -IT_0156 + -IT_0163 + -IT_0166;
    const ccomplex_t IT_0707 = IT_0705 + IT_0706;
    const ccomplex_t IT_0708 = IT_0147*IT_0707;
    const ccomplex_t IT_0709 = IT_0179 + IT_0182 + IT_0190 + IT_0192;
    const ccomplex_t IT_0710 = -IT_0185 + -IT_0187 + -IT_0195;
    const ccomplex_t IT_0711 = IT_0709 + IT_0710;
    const ccomplex_t IT_0712 = IT_0175*IT_0711;
    const ccomplex_t IT_0713 = 2*IT_0208;
    const ccomplex_t IT_0714 = IT_0202 + IT_0203 + IT_0205 + IT_0207 + IT_0210
       + IT_0211 + IT_0212;
    const ccomplex_t IT_0715 = IT_0713 + IT_0714;
    const ccomplex_t IT_0716 = IT_0200*IT_0715;
    const ccomplex_t IT_0717 = 2*IT_0045;
    const ccomplex_t IT_0718 = IT_0039 + IT_0041 + IT_0043 + IT_0047 + IT_0049
       + IT_0050 + IT_0051;
    const ccomplex_t IT_0719 = IT_0717 + IT_0718;
    const ccomplex_t IT_0720 = IT_0222*IT_0719;
    const ccomplex_t IT_0721 = 2*IT_0085;
    const ccomplex_t IT_0722 = IT_0076 + IT_0078 + IT_0080 + IT_0083 + IT_0084
       + IT_0086 + IT_0088;
    const ccomplex_t IT_0723 = IT_0721 + IT_0722;
    const ccomplex_t IT_0724 = IT_0254*IT_0723;
    const ccomplex_t IT_0725 = IT_0278 + IT_0280 + IT_0282 + IT_0284 + IT_0285
       + IT_0287 + IT_0289;
    const ccomplex_t IT_0726 = 2*IT_0286;
    const ccomplex_t IT_0727 = IT_0725 + IT_0726;
    const ccomplex_t IT_0728 = IT_0277*IT_0727;
    const ccomplex_t IT_0729 = 2*IT_0313;
    const ccomplex_t IT_0730 = IT_0308 + IT_0309 + IT_0311 + IT_0315 + IT_0317
       + IT_0318 + IT_0319;
    const ccomplex_t IT_0731 = IT_0729 + IT_0730;
    const ccomplex_t IT_0732 = IT_0307*IT_0731;
    const ccomplex_t IT_0733 = IT_0029 + IT_0103 + IT_0220 + IT_0250 + IT_0273
       + IT_0303 + IT_0700 + IT_0704 + IT_0708 + IT_0712 + IT_0716 + IT_0720 +
       IT_0724 + IT_0728 + IT_0732;
    const ccomplex_t IT_0734 = IT_0039 + IT_0041 + IT_0043 + IT_0047 + IT_0049
       + IT_0050 + IT_0051 + IT_0342;
    const ccomplex_t IT_0735 = IT_0717 + IT_0734;
    const ccomplex_t IT_0736 = IT_0341*IT_0735;
    const ccomplex_t IT_0737 = IT_0370 + IT_0374 + IT_0377;
    const ccomplex_t IT_0738 = -IT_0372 + -IT_0379 + -IT_0382;
    const ccomplex_t IT_0739 = IT_0737 + IT_0738;
    const ccomplex_t IT_0740 = IT_0367*IT_0739;
    const ccomplex_t IT_0741 = IT_0400 + IT_0404 + IT_0407;
    const ccomplex_t IT_0742 = -IT_0402 + -IT_0409 + -IT_0412;
    const ccomplex_t IT_0743 = IT_0741 + IT_0742;
    const ccomplex_t IT_0744 = IT_0397*IT_0743;
    const ccomplex_t IT_0745 = -IT_0432 + -IT_0435 + -IT_0442;
    const ccomplex_t IT_0746 = IT_0430 + IT_0436 + IT_0439;
    const ccomplex_t IT_0747 = IT_0745 + IT_0746;
    const ccomplex_t IT_0748 = IT_0427*IT_0747;
    const ccomplex_t IT_0749 = IT_0447*IT_0739;
    const ccomplex_t IT_0750 = IT_0453*IT_0743;
    const ccomplex_t IT_0751 = IT_0459*IT_0747;
    const ccomplex_t IT_0752 = IT_0225*IT_0715;
    const ccomplex_t IT_0753 = IT_0059*IT_0719;
    const ccomplex_t IT_0754 = IT_0076 + IT_0078 + IT_0080 + IT_0083 + IT_0084
       + IT_0086 + IT_0088 + IT_0466;
    const ccomplex_t IT_0755 = IT_0721 + IT_0754;
    const ccomplex_t IT_0756 = IT_0465*IT_0755;
    const ccomplex_t IT_0757 = IT_0257*IT_0699;
    const ccomplex_t IT_0758 = IT_0073*IT_0723;
    const ccomplex_t IT_0759 = IT_0308 + IT_0309 + IT_0311 + IT_0315 + IT_0317
       + IT_0318 + IT_0319 + IT_0472;
    const ccomplex_t IT_0760 = IT_0729 + IT_0759;
    const ccomplex_t IT_0761 = IT_0471*IT_0760;
    const ccomplex_t IT_0762 = IT_0326*IT_0727;
    const ccomplex_t IT_0763 = IT_0334*IT_0731;
    const ccomplex_t IT_0764 = IT_0124 + IT_0130 + IT_0490;
    const ccomplex_t IT_0765 = -IT_0120 + -IT_0127 + -IT_0489;
    const ccomplex_t IT_0766 = IT_0764 + IT_0765;
    const ccomplex_t IT_0767 = IT_0486*IT_0766;
    const ccomplex_t IT_0768 = IT_0150 + IT_0163 + IT_0506;
    const ccomplex_t IT_0769 = -IT_0153 + -IT_0156 + -IT_0508;
    const ccomplex_t IT_0770 = IT_0768 + IT_0769;
    const ccomplex_t IT_0771 = IT_0504*IT_0770;
    const ccomplex_t IT_0772 = IT_0185 + IT_0190 + IT_0526;
    const ccomplex_t IT_0773 = -IT_0182 + -IT_0187 + -IT_0525;
    const ccomplex_t IT_0774 = IT_0772 + IT_0773;
    const ccomplex_t IT_0775 = IT_0522*IT_0774;
    const ccomplex_t IT_0776 = IT_0531*IT_0766;
    const ccomplex_t IT_0777 = IT_0537*IT_0770;
    const ccomplex_t IT_0778 = IT_0543*IT_0774;
    const ccomplex_t IT_0779 = -IT_0036 + -IT_0067 + -IT_0110 + -IT_0261 + 
      -IT_0330 + -IT_0338 + -IT_0349 + -IT_0352 + -IT_0355 + 0.5*IT_0365 + 0.5
      *IT_0395 + 0.5*IT_0425 + (-0.5)*IT_0451 + (-0.5)*IT_0457 + (-0.5)*IT_0463 
      + 0.5*IT_0484 + 0.5*IT_0502 + 0.5*IT_0520 + (-0.5)*IT_0535 + (-0.5)
      *IT_0541 + (-0.5)*IT_0547 + 2*IT_0736 + (-0.5)*IT_0740 + (-0.5)*IT_0744 + 
      (-0.5)*IT_0748 + 0.5*IT_0749 + 0.5*IT_0750 + 0.5*IT_0751 + -IT_0752 + 
      -IT_0753 + 2*IT_0756 + -IT_0757 + -IT_0758 + 2*IT_0761 + -IT_0762 + 
      -IT_0763 + (-0.5)*IT_0767 + (-0.5)*IT_0771 + (-0.5)*IT_0775 + 0.5*IT_0776 
      + 0.5*IT_0777 + 0.5*IT_0778;
    const ccomplex_t IT_0780 = IT_0733 + IT_0779;
    const ccomplex_t IT_0781 = IT_0005*IT_0780;
    const ccomplex_t IT_0782 = 1.4142135623731*IT_0781;
    const ccomplex_t IT_0783 = (-0.25)*IT_0782;
    const ccomplex_t IT_0784 = IT_0203 + IT_0205 + IT_0208 + IT_0211 + IT_0212;
    const ccomplex_t IT_0785 = IT_0200*IT_0784;
    const ccomplex_t IT_0786 = IT_0229 + IT_0231 + IT_0235 + IT_0237 + IT_0238;
    const ccomplex_t IT_0787 = IT_0228*IT_0786;
    const ccomplex_t IT_0788 = IT_0076 + IT_0080 + IT_0084 + IT_0085 + IT_0086;
    const ccomplex_t IT_0789 = IT_0254*IT_0788;
    const ccomplex_t IT_0790 = IT_0309 + IT_0313 + IT_0315 + IT_0317 + IT_0319;
    const ccomplex_t IT_0791 = IT_0307*IT_0790;
    const ccomplex_t IT_0792 = IT_0124 + IT_0490 + IT_0575 + IT_0576;
    const ccomplex_t IT_0793 = (-2)*IT_0127 + -IT_0489 + -IT_0573 + -IT_0579;
    const ccomplex_t IT_0794 = IT_0792 + IT_0793;
    const ccomplex_t IT_0795 = IT_0116*IT_0794;
    const ccomplex_t IT_0796 = IT_0150 + IT_0506 + IT_0588 + IT_0589;
    const ccomplex_t IT_0797 = (-2)*IT_0156 + -IT_0508 + -IT_0583 + -IT_0587;
    const ccomplex_t IT_0798 = IT_0796 + IT_0797;
    const ccomplex_t IT_0799 = IT_0147*IT_0798;
    const ccomplex_t IT_0800 = IT_0190 + IT_0526 + IT_0595 + IT_0599;
    const ccomplex_t IT_0801 = (-2)*IT_0187 + -IT_0525 + -IT_0593 + -IT_0598;
    const ccomplex_t IT_0802 = IT_0800 + IT_0801;
    const ccomplex_t IT_0803 = IT_0175*IT_0802;
    const ccomplex_t IT_0804 = IT_0041 + IT_0045 + IT_0047 + IT_0049 + IT_0051;
    const ccomplex_t IT_0805 = IT_0222*IT_0804;
    const ccomplex_t IT_0806 = IT_0278 + IT_0280 + IT_0284 + IT_0285 + IT_0286;
    const ccomplex_t IT_0807 = IT_0277*IT_0806;
    const ccomplex_t IT_0808 = IT_0553 + IT_0560 + IT_0571 + IT_0572 + IT_0607
       + IT_0625 + IT_0785 + IT_0787 + IT_0789 + IT_0791 + IT_0795 + IT_0799 +
       IT_0803 + IT_0805 + IT_0807;
    const ccomplex_t IT_0809 = IT_0257*IT_0786;
    const ccomplex_t IT_0810 = IT_0073*IT_0788;
    const ccomplex_t IT_0811 = IT_0308 + IT_0309 + IT_0313 + IT_0315 + IT_0317
       + IT_0319;
    const ccomplex_t IT_0812 = IT_0471*IT_0811;
    const ccomplex_t IT_0813 = IT_0334*IT_0790;
    const ccomplex_t IT_0814 = IT_0225*IT_0784;
    const ccomplex_t IT_0815 = IT_0059*IT_0804;
    const ccomplex_t IT_0816 = IT_0124 + IT_0579;
    const ccomplex_t IT_0817 = -IT_0127 + -IT_0576;
    const ccomplex_t IT_0818 = IT_0816 + IT_0817;
    const ccomplex_t IT_0819 = IT_0486*IT_0818;
    const ccomplex_t IT_0820 = IT_0370 + IT_0644;
    const ccomplex_t IT_0821 = -IT_0372 + -IT_0645;
    const ccomplex_t IT_0822 = IT_0820 + IT_0821;
    const ccomplex_t IT_0823 = IT_0367*IT_0822;
    const ccomplex_t IT_0824 = IT_0400 + IT_0652;
    const ccomplex_t IT_0825 = -IT_0402 + -IT_0653;
    const ccomplex_t IT_0826 = IT_0824 + IT_0825;
    const ccomplex_t IT_0827 = IT_0397*IT_0826;
    const ccomplex_t IT_0828 = IT_0430 + IT_0658;
    const ccomplex_t IT_0829 = -IT_0432 + -IT_0659;
    const ccomplex_t IT_0830 = IT_0828 + IT_0829;
    const ccomplex_t IT_0831 = IT_0427*IT_0830;
    const ccomplex_t IT_0832 = IT_0447*IT_0822;
    const ccomplex_t IT_0833 = IT_0453*IT_0826;
    const ccomplex_t IT_0834 = IT_0459*IT_0830;
    const ccomplex_t IT_0835 = IT_0039 + IT_0041 + IT_0045 + IT_0047 + IT_0049
       + IT_0051;
    const ccomplex_t IT_0836 = IT_0341*IT_0835;
    const ccomplex_t IT_0837 = IT_0076 + IT_0080 + IT_0083 + IT_0084 + IT_0085
       + IT_0086;
    const ccomplex_t IT_0838 = IT_0465*IT_0837;
    const ccomplex_t IT_0839 = IT_0326*IT_0806;
    const ccomplex_t IT_0840 = IT_0150 + IT_0587;
    const ccomplex_t IT_0841 = -IT_0156 + -IT_0589;
    const ccomplex_t IT_0842 = IT_0840 + IT_0841;
    const ccomplex_t IT_0843 = IT_0504*IT_0842;
    const ccomplex_t IT_0844 = IT_0190 + IT_0598;
    const ccomplex_t IT_0845 = -IT_0187 + -IT_0595;
    const ccomplex_t IT_0846 = IT_0844 + IT_0845;
    const ccomplex_t IT_0847 = IT_0522*IT_0846;
    const ccomplex_t IT_0848 = IT_0531*IT_0818;
    const ccomplex_t IT_0849 = IT_0537*IT_0842;
    const ccomplex_t IT_0850 = IT_0543*IT_0846;
    const ccomplex_t IT_0851 = -IT_0554 + -IT_0559 + -IT_0563 + -IT_0566 + 
      -IT_0569 + -IT_0570 + -IT_0620 + -IT_0631 + -IT_0633 + 0.5*IT_0636 + (-0.5
      )*IT_0637 + 0.5*IT_0639 + 0.5*IT_0641 + 0.5*IT_0649 + (-0.5)*IT_0663 + (
      -0.5)*IT_0665 + 0.5*IT_0679 + 0.5*IT_0683 + (-0.5)*IT_0687 + (-0.5)
      *IT_0689 + (-0.5)*IT_0691 + -IT_0809 + -IT_0810 + 2*IT_0812 + -IT_0813 + 
      -IT_0814 + -IT_0815 + (-0.5)*IT_0819 + (-0.5)*IT_0823 + (-0.5)*IT_0827 + (
      -0.5)*IT_0831 + 0.5*IT_0832 + 0.5*IT_0833 + 0.5*IT_0834 + 2*IT_0836 + 2
      *IT_0838 + -IT_0839 + (-0.5)*IT_0843 + (-0.5)*IT_0847 + 0.5*IT_0848 + 0.5
      *IT_0849 + 0.5*IT_0850;
    const ccomplex_t IT_0852 = IT_0808 + IT_0851;
    const ccomplex_t IT_0853 = IT_0005*IT_0852;
    const ccomplex_t IT_0854 = 1.4142135623731*IT_0853;
    const ccomplex_t IT_0855 = 0.25*IT_0854;
    return create_ccomplex_return((0 + _Complex_I*0.25)*IT_0552 + (0 +
       _Complex_I*0.25)*IT_0696 + (0 + _Complex_I*0.25)*IT_0783 + (0 +
       _Complex_I*0.25)*IT_0855);
}

