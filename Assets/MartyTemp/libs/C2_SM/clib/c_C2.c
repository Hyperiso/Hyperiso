#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "marty/core/looptools_interface.h"
#include "c_C2.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_C2(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t G_F = param->G_F;
    const creal_t M_W = param->M_W;
    const creal_t M_Z = param->M_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t s_13 = param->s_13;
    const creal_t s_24 = param->s_24;
    const creal_t Finite = param->Finite;
    const creal_t theta_W = param->theta_W;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = pow(G_F, -1);
    const ccomplex_t IT_0001 = pow(V_cb, -1);
    const ccomplex_t IT_0002 = cpow(conj(V_cs), -1);
    const ccomplex_t IT_0003 = pow(M_W, 2);
    const ccomplex_t IT_0004 = pow(m_b, 2);
    const ccomplex_t IT_0005 = pow(m_c, 2);
    const ccomplex_t IT_0006 = cpow(2*s_13 + IT_0003 + -IT_0004 + -IT_0005 + 
      -reg_prop, -1);
    const ccomplex_t IT_0007 = pow(1.4142135623731, 0.5);
    const ccomplex_t IT_0008 = pow(G_F, 0.5);
    const ccomplex_t IT_0009 = (0 + _Complex_I*2.82842712474619)*M_W*conj(V_cs
      )*IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*2.82842712474619)*M_W*V_cb
      *IT_0007*IT_0008;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = pow(M_Z, 2);
    const ccomplex_t IT_0017 = pow(m_s, 2);
    const ccomplex_t IT_0018 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0016, IT_0005, IT_0017);
    const ccomplex_t IT_0019 = (0 + _Complex_I*2.82842712474619)*m_s*conj(V_cs
      )*IT_0007*IT_0008;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = 2*m_c*IT_0007*IT_0008;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = 2*m_s*IT_0007*IT_0008;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*2.82842712474619)*m_c*V_cb
      *IT_0007*IT_0008;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0020*IT_0022*IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*0.101321183642338)*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*0.101321183642338)*m_b*m_c;
    const ccomplex_t IT_0030 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0016, IT_0004, IT_0005);
    const ccomplex_t IT_0031 = cos(theta_W);
    const ccomplex_t IT_0032 = cpow(IT_0031, -1);
    const ccomplex_t IT_0033 = sin(theta_W);
    const ccomplex_t IT_0034 = cpow(IT_0033, 2);
    const ccomplex_t IT_0035 = (0 + _Complex_I*2)*M_W*IT_0007*IT_0008*IT_0032
      *IT_0034;
    const ccomplex_t IT_0036 = (-0.666666666666667)*IT_0035;
    const ccomplex_t IT_0037 = 0.333333333333333*IT_0035;
    const ccomplex_t IT_0038 = IT_0020*IT_0026*IT_0030*IT_0036*IT_0037;
    const ccomplex_t IT_0039 = IT_0029*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0037, 2);
    const ccomplex_t IT_0041 = (0 + _Complex_I*0.101321183642338)*m_b*m_s
      *IT_0040;
    const ccomplex_t IT_0042 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0016, IT_0004, IT_0017);
    const ccomplex_t IT_0043 = (0 + _Complex_I*2.82842712474619)*m_c*conj(V_cs
      )*IT_0007*IT_0008;
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = IT_0026*IT_0042*IT_0044;
    const ccomplex_t IT_0046 = IT_0041*IT_0045;
    const ccomplex_t IT_0047 = IT_0000*IT_0001*IT_0002;
    const ccomplex_t IT_0048 = IT_0039 + IT_0046;
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = m_c*m_s;
    const ccomplex_t IT_0052 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0016, IT_0005, IT_0017);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = IT_0010*IT_0012*IT_0022*IT_0024;
    const ccomplex_t IT_0055 = (0 + _Complex_I*0.101321183642338)*IT_0054;
    const ccomplex_t IT_0056 = IT_0022*IT_0024*IT_0026*IT_0044;
    const ccomplex_t IT_0057 = (0 + _Complex_I*0.101321183642338)*IT_0056;
    const ccomplex_t IT_0058 = IT_0032*IT_0033;
    const ccomplex_t IT_0059 = IT_0007*IT_0008*IT_0033;
    const ccomplex_t IT_0060 = 4*M_W*IT_0058*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = cpow(IT_0033, -1);
    const ccomplex_t IT_0063 = IT_0031*IT_0062;
    const ccomplex_t IT_0064 = 4*M_W*IT_0059*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0061 + 3*IT_0065);
    const ccomplex_t IT_0067 = (-0.166666666666667)*IT_0066;
    const ccomplex_t IT_0068 = cpow(IT_0067, 2);
    const ccomplex_t IT_0069 = (0 + _Complex_I*0.101321183642338)*IT_0068;
    const ccomplex_t IT_0070 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0016, IT_0004, IT_0017);
    const ccomplex_t IT_0071 = IT_0010*IT_0012*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0061 + (-3)*IT_0065);
    const ccomplex_t IT_0074 = (-0.166666666666667)*IT_0073;
    const ccomplex_t IT_0075 = cpow(IT_0074, 2);
    const ccomplex_t IT_0076 = (0 + _Complex_I*0.101321183642338)*IT_0075;
    const ccomplex_t IT_0077 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0016, IT_0005, IT_0005);
    const ccomplex_t IT_0078 = IT_0010*IT_0012*IT_0077;
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0004, IT_0017);
    const ccomplex_t IT_0081 = (0 + _Complex_I*2)*M_W*IT_0007*IT_0008*IT_0033;
    const ccomplex_t IT_0082 = (-0.333333333333333)*IT_0081;
    const ccomplex_t IT_0083 = cpow(IT_0082, 2);
    const ccomplex_t IT_0084 = (0 + _Complex_I*0.101321183642338)*IT_0083;
    const ccomplex_t IT_0085 = IT_0013*IT_0084;
    const ccomplex_t IT_0086 = IT_0080*IT_0085;
    const ccomplex_t IT_0087 = IT_0072 + IT_0079 + IT_0086;
    const ccomplex_t IT_0088 = IT_0010*IT_0012*IT_0018*IT_0067*IT_0074;
    const ccomplex_t IT_0089 = (0 + _Complex_I*0.101321183642338)*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*0.101321183642338)*m_c*m_s;
    const ccomplex_t IT_0091 = IT_0020*IT_0026*IT_0052*IT_0067*IT_0074;
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0016, IT_0004, IT_0005);
    const ccomplex_t IT_0094 = IT_0010*IT_0012*IT_0067*IT_0074*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*0.101321183642338)*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*2.82842712474619)*m_b*V_cb
      *IT_0007*IT_0008;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = IT_0030*IT_0044*IT_0067*IT_0074*IT_0097;
    const ccomplex_t IT_0099 = IT_0029*IT_0098;
    const ccomplex_t IT_0100 = IT_0026*IT_0030*IT_0037*IT_0044*IT_0074;
    const ccomplex_t IT_0101 = IT_0029*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*0.101321183642338)*IT_0005;
    const ccomplex_t IT_0103 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0016, IT_0005, IT_0005);
    const ccomplex_t IT_0104 = IT_0020*IT_0026*IT_0036*IT_0074*IT_0103;
    const ccomplex_t IT_0105 = IT_0102*IT_0104;
    const ccomplex_t IT_0106 = IT_0026*IT_0037*IT_0044*IT_0052*IT_0074;
    const ccomplex_t IT_0107 = IT_0090*IT_0106;
    const ccomplex_t IT_0108 = IT_0036*IT_0044*IT_0074*IT_0097*IT_0103;
    const ccomplex_t IT_0109 = IT_0102*IT_0108;
    const ccomplex_t IT_0110 = IT_0020*IT_0042*IT_0097;
    const ccomplex_t IT_0111 = (0 + _Complex_I*0.101321183642338)*m_b*m_s
      *IT_0068;
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*0.101321183642338)*m_b*m_s;
    const ccomplex_t IT_0114 = IT_0020*IT_0026*IT_0037*IT_0042*IT_0067;
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = IT_0020*IT_0036*IT_0052*IT_0067*IT_0097;
    const ccomplex_t IT_0117 = IT_0090*IT_0116;
    const ccomplex_t IT_0118 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0005, IT_0005);
    const ccomplex_t IT_0119 = IT_0005*IT_0118;
    const ccomplex_t IT_0120 = IT_0044*IT_0097;
    const ccomplex_t IT_0121 = 0.666666666666667*IT_0081;
    const ccomplex_t IT_0122 = cpow(IT_0121, 2);
    const ccomplex_t IT_0123 = (0 + _Complex_I*0.101321183642338)*IT_0122;
    const ccomplex_t IT_0124 = IT_0120*IT_0123;
    const ccomplex_t IT_0125 = IT_0119*IT_0124;
    const ccomplex_t IT_0126 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0005, IT_0017);
    const ccomplex_t IT_0127 = IT_0010*IT_0012*IT_0082*IT_0121;
    const ccomplex_t IT_0128 = (0 + _Complex_I*0.101321183642338)*IT_0127;
    const ccomplex_t IT_0129 = IT_0126*IT_0128;
    const ccomplex_t IT_0130 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0005, IT_0017);
    const ccomplex_t IT_0131 = IT_0051*IT_0130;
    const ccomplex_t IT_0132 = IT_0020*IT_0082*IT_0097*IT_0121;
    const ccomplex_t IT_0133 = (0 + _Complex_I*0.101321183642338)*IT_0132;
    const ccomplex_t IT_0134 = IT_0131*IT_0133;
    const ccomplex_t IT_0135 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0004, IT_0005);
    const ccomplex_t IT_0136 = m_b*m_c;
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = IT_0133*IT_0137;
    const ccomplex_t IT_0139 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0004, IT_0005);
    const ccomplex_t IT_0140 = IT_0128*IT_0139;
    const ccomplex_t IT_0141 = IT_0020*IT_0022*IT_0024*IT_0097;
    const ccomplex_t IT_0142 = (0 + _Complex_I*0.101321183642338)*IT_0141;
    const ccomplex_t IT_0143 = IT_0018*IT_0142;
    const ccomplex_t IT_0144 = IT_0026*IT_0044*IT_0082*IT_0121;
    const ccomplex_t IT_0145 = (0 + _Complex_I*0.101321183642338)*IT_0144;
    const ccomplex_t IT_0146 = IT_0137*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*0.101321183642338)*IT_0005
      *IT_0075;
    const ccomplex_t IT_0148 = IT_0026*IT_0044*IT_0103;
    const ccomplex_t IT_0149 = IT_0147*IT_0148;
    const ccomplex_t IT_0150 = IT_0131*IT_0145;
    const ccomplex_t IT_0151 = IT_0020*IT_0026*IT_0082*IT_0121;
    const ccomplex_t IT_0152 = (0 + _Complex_I*0.101321183642338)*IT_0151;
    const ccomplex_t IT_0153 = IT_0131*IT_0152;
    const ccomplex_t IT_0154 = IT_0044*IT_0082*IT_0097*IT_0121;
    const ccomplex_t IT_0155 = (0 + _Complex_I*0.101321183642338)*IT_0154;
    const ccomplex_t IT_0156 = IT_0131*IT_0155;
    const ccomplex_t IT_0157 = IT_0137*IT_0152;
    const ccomplex_t IT_0158 = IT_0137*IT_0155;
    const ccomplex_t IT_0159 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0004, IT_0017);
    const ccomplex_t IT_0160 = m_b*m_s;
    const ccomplex_t IT_0161 = IT_0159*IT_0160;
    const ccomplex_t IT_0162 = IT_0020*IT_0026;
    const ccomplex_t IT_0163 = IT_0084*IT_0162;
    const ccomplex_t IT_0164 = IT_0161*IT_0163;
    const ccomplex_t IT_0165 = IT_0026*IT_0044;
    const ccomplex_t IT_0166 = IT_0084*IT_0165;
    const ccomplex_t IT_0167 = IT_0161*IT_0166;
    const ccomplex_t IT_0168 = IT_0020*IT_0097;
    const ccomplex_t IT_0169 = IT_0084*IT_0168;
    const ccomplex_t IT_0170 = IT_0161*IT_0169;
    const ccomplex_t IT_0171 = IT_0084*IT_0120;
    const ccomplex_t IT_0172 = IT_0161*IT_0171;
    const ccomplex_t IT_0173 = pow(m_h, 2);
    const ccomplex_t IT_0174 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0004, IT_0005, IT_0173);
    const ccomplex_t IT_0175 = IT_0136*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*2)*m_c*IT_0007*IT_0008;
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*2)*m_b*IT_0007*IT_0008;
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = IT_0010*IT_0012*IT_0177*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*0.101321183642338)*IT_0180;
    const ccomplex_t IT_0182 = IT_0175*IT_0181;
    const ccomplex_t IT_0183 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0004, IT_0173, IT_0017);
    const ccomplex_t IT_0184 = IT_0160*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*2)*m_s*IT_0007*IT_0008;
    const ccomplex_t IT_0186 = (-0.5)*IT_0185;
    const ccomplex_t IT_0187 = IT_0010*IT_0012*IT_0179*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*0.101321183642338)*IT_0187;
    const ccomplex_t IT_0189 = IT_0184*IT_0188;
    const ccomplex_t IT_0190 = IT_0030*IT_0136;
    const ccomplex_t IT_0191 = 2*m_b*IT_0007*IT_0008;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = IT_0010*IT_0012*IT_0022*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*0.101321183642338)*IT_0193;
    const ccomplex_t IT_0195 = IT_0190*IT_0194;
    const ccomplex_t IT_0196 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0005, IT_0005, IT_0173);
    const ccomplex_t IT_0197 = IT_0005*IT_0196;
    const ccomplex_t IT_0198 = cpow(IT_0177, 2);
    const ccomplex_t IT_0199 = (0 + _Complex_I*0.101321183642338)*IT_0198;
    const ccomplex_t IT_0200 = IT_0013*IT_0199;
    const ccomplex_t IT_0201 = IT_0197*IT_0200;
    const ccomplex_t IT_0202 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0005, IT_0173, IT_0017);
    const ccomplex_t IT_0203 = IT_0051*IT_0202;
    const ccomplex_t IT_0204 = IT_0010*IT_0012*IT_0177*IT_0186;
    const ccomplex_t IT_0205 = (0 + _Complex_I*0.101321183642338)*IT_0204;
    const ccomplex_t IT_0206 = IT_0203*IT_0205;
    const ccomplex_t IT_0207 = IT_0005*IT_0103;
    const ccomplex_t IT_0208 = cpow(IT_0022, 2);
    const ccomplex_t IT_0209 = (0 + _Complex_I*0.101321183642338)*IT_0208;
    const ccomplex_t IT_0210 = IT_0013*IT_0209;
    const ccomplex_t IT_0211 = IT_0207*IT_0210;
    const ccomplex_t IT_0212 = IT_0018*IT_0057;
    const ccomplex_t IT_0213 = IT_0022*IT_0024*IT_0044*IT_0097;
    const ccomplex_t IT_0214 = (0 + _Complex_I*0.101321183642338)*IT_0213;
    const ccomplex_t IT_0215 = IT_0018*IT_0214;
    const ccomplex_t IT_0216 = (-4)*IT_0089 + (-0.25)*IT_0092 + (-4)*IT_0095 +
       (-0.25)*IT_0099 + (-0.25)*IT_0101 + (-0.25)*IT_0105 + (-0.25)*IT_0107 + (
      -0.25)*IT_0109 + (-0.25)*IT_0112 + (-0.25)*IT_0115 + (-0.25)*IT_0117 + (
      -0.25)*IT_0125 + (-4)*IT_0129 + (-0.25)*IT_0134 + (-0.25)*IT_0138 + (-4)
      *IT_0140 + 0.25*IT_0143 + (-0.25)*IT_0146 + (-0.25)*IT_0149 + (-0.25)
      *IT_0150 + (-0.25)*IT_0153 + (-0.25)*IT_0156 + (-0.25)*IT_0157 + (-0.25)
      *IT_0158 + (-0.25)*IT_0164 + (-0.25)*IT_0167 + (-0.25)*IT_0170 + (-0.25)
      *IT_0172 + (-0.25)*IT_0182 + (-0.25)*IT_0189 + (-0.25)*IT_0195 + (-0.25)
      *IT_0201 + (-0.25)*IT_0206 + 0.25*IT_0211 + 0.25*IT_0212 + (-0.25)*IT_0215;
    const ccomplex_t IT_0217 = IT_0087 + IT_0216;
    const ccomplex_t IT_0218 = IT_0047*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = 0.015625*IT_0219;
    const ccomplex_t IT_0221 = IT_0020*IT_0030*IT_0036*IT_0067*IT_0097;
    const ccomplex_t IT_0222 = IT_0029*IT_0221;
    const ccomplex_t IT_0223 = IT_0037*IT_0042*IT_0044*IT_0067*IT_0097;
    const ccomplex_t IT_0224 = IT_0113*IT_0223;
    const ccomplex_t IT_0225 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0004, IT_0005, IT_0173);
    const ccomplex_t IT_0226 = IT_0020*IT_0097*IT_0177*IT_0179;
    const ccomplex_t IT_0227 = (0 + _Complex_I*0.101321183642338)*IT_0226;
    const ccomplex_t IT_0228 = IT_0225*IT_0227;
    const ccomplex_t IT_0229 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0005, IT_0173, IT_0017);
    const ccomplex_t IT_0230 = IT_0026*IT_0044*IT_0177*IT_0186;
    const ccomplex_t IT_0231 = (0 + _Complex_I*0.101321183642338)*IT_0230;
    const ccomplex_t IT_0232 = IT_0229*IT_0231;
    const ccomplex_t IT_0233 = IT_0168*IT_0209;
    const ccomplex_t IT_0234 = IT_0077*IT_0233;
    const ccomplex_t IT_0235 = IT_0123*IT_0168;
    const ccomplex_t IT_0236 = IT_0119*IT_0235;
    const ccomplex_t IT_0237 = IT_0123*IT_0162;
    const ccomplex_t IT_0238 = IT_0119*IT_0237;
    const ccomplex_t IT_0239 = IT_0123*IT_0165;
    const ccomplex_t IT_0240 = IT_0119*IT_0239;
    const ccomplex_t IT_0241 = IT_0020*IT_0026*IT_0177*IT_0179;
    const ccomplex_t IT_0242 = (0 + _Complex_I*0.101321183642338)*IT_0241;
    const ccomplex_t IT_0243 = IT_0225*IT_0242;
    const ccomplex_t IT_0244 = IT_0026*IT_0044*IT_0177*IT_0179;
    const ccomplex_t IT_0245 = (0 + _Complex_I*0.101321183642338)*IT_0244;
    const ccomplex_t IT_0246 = IT_0225*IT_0245;
    const ccomplex_t IT_0247 = IT_0044*IT_0097*IT_0177*IT_0179;
    const ccomplex_t IT_0248 = (0 + _Complex_I*0.101321183642338)*IT_0247;
    const ccomplex_t IT_0249 = IT_0225*IT_0248;
    const ccomplex_t IT_0250 = IT_0020*IT_0022*IT_0026*IT_0192;
    const ccomplex_t IT_0251 = (0 + _Complex_I*0.101321183642338)*IT_0250;
    const ccomplex_t IT_0252 = IT_0093*IT_0251;
    const ccomplex_t IT_0253 = IT_0022*IT_0044*IT_0097*IT_0192;
    const ccomplex_t IT_0254 = (0 + _Complex_I*0.101321183642338)*IT_0253;
    const ccomplex_t IT_0255 = IT_0093*IT_0254;
    const ccomplex_t IT_0256 = IT_0024*IT_0026*IT_0044*IT_0192;
    const ccomplex_t IT_0257 = (0 + _Complex_I*0.101321183642338)*IT_0256;
    const ccomplex_t IT_0258 = IT_0070*IT_0257;
    const ccomplex_t IT_0259 = IT_0020*IT_0024*IT_0097*IT_0192;
    const ccomplex_t IT_0260 = (0 + _Complex_I*0.101321183642338)*IT_0259;
    const ccomplex_t IT_0261 = IT_0070*IT_0260;
    const ccomplex_t IT_0262 = IT_0020*IT_0026*IT_0177*IT_0186;
    const ccomplex_t IT_0263 = (0 + _Complex_I*0.101321183642338)*IT_0262;
    const ccomplex_t IT_0264 = IT_0229*IT_0263;
    const ccomplex_t IT_0265 = IT_0020*IT_0097*IT_0177*IT_0186;
    const ccomplex_t IT_0266 = (0 + _Complex_I*0.101321183642338)*IT_0265;
    const ccomplex_t IT_0267 = IT_0229*IT_0266;
    const ccomplex_t IT_0268 = IT_0044*IT_0097*IT_0177*IT_0186;
    const ccomplex_t IT_0269 = (0 + _Complex_I*0.101321183642338)*IT_0268;
    const ccomplex_t IT_0270 = IT_0229*IT_0269;
    const ccomplex_t IT_0271 = IT_0165*IT_0209;
    const ccomplex_t IT_0272 = IT_0077*IT_0271;
    const ccomplex_t IT_0273 = IT_0222 + IT_0224 + IT_0228 + IT_0232 + IT_0234
       + IT_0236 + IT_0238 + IT_0240 + IT_0243 + IT_0246 + IT_0249 + IT_0252 +
       IT_0255 + IT_0258 + IT_0261 + IT_0264 + IT_0267 + IT_0270 + IT_0272;
    const ccomplex_t IT_0274 = IT_0042*IT_0160;
    const ccomplex_t IT_0275 = IT_0010*IT_0012*IT_0024*IT_0192;
    const ccomplex_t IT_0276 = (0 + _Complex_I*0.101321183642338)*IT_0275;
    const ccomplex_t IT_0277 = IT_0274*IT_0276;
    const ccomplex_t IT_0278 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0005, IT_0005, IT_0173);
    const ccomplex_t IT_0279 = IT_0165*IT_0199;
    const ccomplex_t IT_0280 = IT_0278*IT_0279;
    const ccomplex_t IT_0281 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0004, IT_0173, IT_0017);
    const ccomplex_t IT_0282 = IT_0020*IT_0026*IT_0179*IT_0186;
    const ccomplex_t IT_0283 = (0 + _Complex_I*0.101321183642338)*IT_0282;
    const ccomplex_t IT_0284 = IT_0281*IT_0283;
    const ccomplex_t IT_0285 = IT_0020*IT_0097*IT_0179*IT_0186;
    const ccomplex_t IT_0286 = (0 + _Complex_I*0.101321183642338)*IT_0285;
    const ccomplex_t IT_0287 = IT_0281*IT_0286;
    const ccomplex_t IT_0288 = IT_0162*IT_0199;
    const ccomplex_t IT_0289 = IT_0278*IT_0288;
    const ccomplex_t IT_0290 = IT_0020*IT_0022*IT_0097*IT_0192;
    const ccomplex_t IT_0291 = (0 + _Complex_I*0.101321183642338)*IT_0290;
    const ccomplex_t IT_0292 = IT_0093*IT_0291;
    const ccomplex_t IT_0293 = IT_0024*IT_0044*IT_0097*IT_0192;
    const ccomplex_t IT_0294 = (0 + _Complex_I*0.101321183642338)*IT_0293;
    const ccomplex_t IT_0295 = IT_0070*IT_0294;
    const ccomplex_t IT_0296 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0003,
       IT_0005, IT_0005);
    const ccomplex_t IT_0297 = IT_0013*IT_0123;
    const ccomplex_t IT_0298 = IT_0296*IT_0297;
    const ccomplex_t IT_0299 = IT_0026*IT_0044*IT_0179*IT_0186;
    const ccomplex_t IT_0300 = (0 + _Complex_I*0.101321183642338)*IT_0299;
    const ccomplex_t IT_0301 = IT_0281*IT_0300;
    const ccomplex_t IT_0302 = IT_0044*IT_0097*IT_0179*IT_0186;
    const ccomplex_t IT_0303 = (0 + _Complex_I*0.101321183642338)*IT_0302;
    const ccomplex_t IT_0304 = IT_0281*IT_0303;
    const ccomplex_t IT_0305 = IT_0022*IT_0026*IT_0044*IT_0192;
    const ccomplex_t IT_0306 = (0 + _Complex_I*0.101321183642338)*IT_0305;
    const ccomplex_t IT_0307 = IT_0093*IT_0306;
    const ccomplex_t IT_0308 = IT_0020*IT_0024*IT_0026*IT_0192;
    const ccomplex_t IT_0309 = (0 + _Complex_I*0.101321183642338)*IT_0308;
    const ccomplex_t IT_0310 = IT_0070*IT_0309;
    const ccomplex_t IT_0311 = IT_0168*IT_0199;
    const ccomplex_t IT_0312 = IT_0278*IT_0311;
    const ccomplex_t IT_0313 = IT_0120*IT_0199;
    const ccomplex_t IT_0314 = IT_0278*IT_0313;
    const ccomplex_t IT_0315 = IT_0162*IT_0209;
    const ccomplex_t IT_0316 = IT_0077*IT_0315;
    const ccomplex_t IT_0317 = IT_0120*IT_0209;
    const ccomplex_t IT_0318 = IT_0077*IT_0317;
    const ccomplex_t IT_0319 = -IT_0277 + -IT_0280 + -IT_0284 + -IT_0287 + 
      -IT_0289 + -IT_0292 + -IT_0295 + (-4)*IT_0298 + -IT_0301 + -IT_0304 + 
      -IT_0307 + -IT_0310 + -IT_0312 + -IT_0314 + -IT_0316 + -IT_0318;
    const ccomplex_t IT_0320 = IT_0273 + IT_0319;
    const ccomplex_t IT_0321 = IT_0047*IT_0320;
    const ccomplex_t IT_0322 = 1.4142135623731*IT_0321;
    const ccomplex_t IT_0323 = (-0.00390625)*IT_0322;
    const ccomplex_t IT_0324 = IT_0072 + IT_0079;
    const ccomplex_t IT_0325 = (-4)*IT_0089 + (-4)*IT_0095 + (-0.25)*IT_0099 +
       0.25*IT_0101 + (-0.25)*IT_0105 + 0.25*IT_0109 + (-0.25)*IT_0112 + (-0.25)
      *IT_0143 + 0.25*IT_0146 + 0.25*IT_0164 + (-0.25)*IT_0172 + (-0.25)*IT_0182
       + (-0.25)*IT_0201 + 0.25*IT_0215 + (-0.25)*IT_0222;
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = IT_0047*IT_0326;
    const ccomplex_t IT_0328 = 1.4142135623731*IT_0327;
    const ccomplex_t IT_0329 = (-0.015625)*IT_0328;
    const ccomplex_t IT_0330 = IT_0053*IT_0055;
    const ccomplex_t IT_0331 = IT_0018*IT_0028;
    const ccomplex_t IT_0332 = IT_0092 + IT_0107 + IT_0138 + IT_0149 + IT_0150
       + IT_0153 + IT_0158 + IT_0170 + IT_0189 + IT_0195 + IT_0206 + IT_0224 +
       IT_0228 + IT_0232 + IT_0238 + IT_0240 + IT_0249 + IT_0255 + IT_0261 +
       IT_0264 + IT_0272 + IT_0284 + IT_0301 + IT_0307 + IT_0310 + IT_0312 +
       IT_0314 + IT_0318 + IT_0330 + IT_0331;
    const ccomplex_t IT_0333 = (-4)*IT_0086 + -IT_0115 + -IT_0117 + -IT_0125 +
       16*IT_0129 + -IT_0134 + 16*IT_0140 + -IT_0156 + -IT_0157 + -IT_0167 + 
      -IT_0211 + -IT_0234 + -IT_0236 + -IT_0243 + -IT_0246 + -IT_0252 + -IT_0258
       + -IT_0267 + -IT_0270 + -IT_0277 + -IT_0280 + -IT_0287 + -IT_0289 + 
      -IT_0292 + -IT_0295 + (-4)*IT_0298 + -IT_0304 + -IT_0316;
    const ccomplex_t IT_0334 = IT_0332 + IT_0333;
    const ccomplex_t IT_0335 = IT_0047*IT_0334;
    const ccomplex_t IT_0336 = 1.4142135623731*IT_0335;
    const ccomplex_t IT_0337 = 0.00390625*IT_0336;
    const ccomplex_t IT_0338 = (-4)*IT_0089 + (-4)*IT_0095 + (-0.25)*IT_0099 +
       (-0.25)*IT_0101 + (-0.25)*IT_0112 + (-4)*IT_0129 + 0.25*IT_0138 + 0.25
      *IT_0150 + 0.25*IT_0156 + 0.25*IT_0157 + 0.25*IT_0167 + 0.25*IT_0172 + (
      -0.25)*IT_0189 + (-0.25)*IT_0206 + 0.25*IT_0215 + 0.25*IT_0222 + 0.25
      *IT_0277 + (-0.25)*IT_0330;
    const ccomplex_t IT_0339 = IT_0087 + IT_0338;
    const ccomplex_t IT_0340 = IT_0047*IT_0339;
    const ccomplex_t IT_0341 = 1.4142135623731*IT_0340;
    const ccomplex_t IT_0342 = (-0.015625)*IT_0341;
    const ccomplex_t IT_0343 = IT_0092 + IT_0109 + IT_0115 + IT_0117 + IT_0125
       + IT_0134 + IT_0146 + IT_0149 + IT_0153 + IT_0158 + IT_0164 + IT_0170 +
       IT_0182 + IT_0195 + IT_0201 + IT_0212 + IT_0240 + IT_0246 + IT_0249 +
       IT_0255 + IT_0261 + IT_0264 + IT_0267 + IT_0272 + IT_0289 + IT_0292 +
       IT_0295 + IT_0301 + IT_0304 + IT_0312 + IT_0316 + IT_0331;
    const ccomplex_t IT_0344 = -IT_0105 + -IT_0107 + 16*IT_0140 + -IT_0143 + 
      -IT_0211 + -IT_0224 + -IT_0228 + -IT_0232 + -IT_0234 + -IT_0236 + -IT_0238
       + -IT_0243 + -IT_0252 + -IT_0258 + -IT_0270 + -IT_0280 + -IT_0284 + 
      -IT_0287 + (-4)*IT_0298 + -IT_0307 + -IT_0310 + -IT_0314 + -IT_0318;
    const ccomplex_t IT_0345 = IT_0343 + IT_0344;
    const ccomplex_t IT_0346 = IT_0047*IT_0345;
    const ccomplex_t IT_0347 = 1.4142135623731*IT_0346;
    const ccomplex_t IT_0348 = 0.00390625*IT_0347;
    const ccomplex_t IT_0349 = (-4)*IT_0089 + (-4)*IT_0095 + 0.25*IT_0101 +
       0.25*IT_0107 + 0.25*IT_0115 + 0.25*IT_0125 + 0.25*IT_0146 + (-0.25)
      *IT_0149 + (-0.25)*IT_0153 + (-0.25)*IT_0156 + (-0.25)*IT_0158 + (-0.25)
      *IT_0167 + (-0.25)*IT_0182 + (-0.25)*IT_0195 + (-0.25)*IT_0206 + (-0.25)
      *IT_0212 + (-0.25)*IT_0215 + 0.25*IT_0222 + 0.25*IT_0232 + (-0.25)*IT_0234
       + 0.25*IT_0246 + (-0.25)*IT_0249 + (-0.25)*IT_0255 + (-0.25)*IT_0261 + (
      -0.25)*IT_0270 + (-0.25)*IT_0284 + 0.25*IT_0287 + (-0.25)*IT_0289 + (-0.25
      )*IT_0307 + (-0.25)*IT_0310 + 0.25*IT_0312 + (-0.25)*IT_0316 + (-0.25)
      *IT_0330;
    const ccomplex_t IT_0350 = IT_0324 + IT_0349;
    const ccomplex_t IT_0351 = IT_0047*IT_0350;
    const ccomplex_t IT_0352 = 1.4142135623731*IT_0351;
    const ccomplex_t IT_0353 = 0.015625*IT_0352;
    const ccomplex_t IT_0354 = IT_0092 + IT_0099 + IT_0112 + IT_0143 + IT_0157
       + IT_0170 + IT_0189 + IT_0201 + IT_0236 + IT_0240 + IT_0243 + IT_0252 +
       IT_0258 + IT_0264 + IT_0272 + IT_0292 + IT_0295 + IT_0304 + IT_0314 +
       IT_0318 + IT_0331;
    const ccomplex_t IT_0355 = (-4)*IT_0086 + -IT_0105 + -IT_0109 + -IT_0117 +
       16*IT_0129 + -IT_0134 + -IT_0138 + 16*IT_0140 + -IT_0150 + -IT_0164 + 
      -IT_0172 + -IT_0211 + -IT_0224 + -IT_0228 + -IT_0238 + -IT_0267 + -IT_0277
       + -IT_0280 + (-4)*IT_0298 + -IT_0301;
    const ccomplex_t IT_0356 = IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = IT_0047*IT_0356;
    const ccomplex_t IT_0358 = 1.4142135623731*IT_0357;
    const ccomplex_t IT_0359 = (-0.00390625)*IT_0358;
    const ccomplex_t IT_0360 = mtylt_C0iC(9, 0, 0, 0, IT_0003, IT_0016,
       IT_0005, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0361 = 2*M_W*IT_0007*IT_0008;
    const ccomplex_t IT_0362 = -IT_0361;
    const ccomplex_t IT_0363 = IT_0360*IT_0362;
    const ccomplex_t IT_0364 = IT_0022*IT_0044*IT_0363;
    const ccomplex_t IT_0365 = 0.101321183642338*IT_0364;
    const ccomplex_t IT_0366 = mtylt_C0iC(9, 0, 0, 0, IT_0003, IT_0016,
       IT_0017, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0367 = IT_0362*IT_0366;
    const ccomplex_t IT_0368 = IT_0024*IT_0044*IT_0367;
    const ccomplex_t IT_0369 = 0.101321183642338*IT_0368;
    const ccomplex_t IT_0370 = mtylt_C0iC(9, 0, 0, 0, IT_0016, IT_0005,
       IT_0017, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0371 = (-2)*IT_0370;
    const ccomplex_t IT_0372 = Finite + IT_0371;
    const ccomplex_t IT_0373 = IT_0010*IT_0067*IT_0074*IT_0372;
    const ccomplex_t IT_0374 = 0.101321183642338*IT_0373;
    const ccomplex_t IT_0375 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0003, IT_0005,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0376 = m_c*IT_0375;
    const ccomplex_t IT_0377 = (0 + _Complex_I*2)*IT_0003*IT_0007*IT_0008
      *IT_0033;
    const ccomplex_t IT_0378 = IT_0020*IT_0121*IT_0376*IT_0377;
    const ccomplex_t IT_0379 = 0.101321183642338*IT_0378;
    const ccomplex_t IT_0380 = IT_0044*IT_0121*IT_0376*IT_0377;
    const ccomplex_t IT_0381 = 0.101321183642338*IT_0380;
    const ccomplex_t IT_0382 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0003, IT_0017,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0383 = m_s*IT_0382;
    const ccomplex_t IT_0384 = IT_0020*IT_0082*IT_0377*IT_0383;
    const ccomplex_t IT_0385 = 0.101321183642338*IT_0384;
    const ccomplex_t IT_0386 = IT_0044*IT_0082*IT_0377*IT_0383;
    const ccomplex_t IT_0387 = 0.101321183642338*IT_0386;
    const ccomplex_t IT_0388 = mtylt_C0iC(0, 0, 0, 0, IT_0003, IT_0005,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0389 = m_c*IT_0388;
    const ccomplex_t IT_0390 = (0 + _Complex_I*2)*IT_0003*IT_0007*IT_0008;
    const ccomplex_t IT_0391 = IT_0010*IT_0177*IT_0389*IT_0390;
    const ccomplex_t IT_0392 = 0.101321183642338*IT_0391;
    const ccomplex_t IT_0393 = mtylt_C0iC(0, 0, 0, 0, IT_0003, IT_0173,
       IT_0017, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0394 = m_s*IT_0393;
    const ccomplex_t IT_0395 = IT_0010*IT_0186*IT_0390*IT_0394;
    const ccomplex_t IT_0396 = 0.101321183642338*IT_0395;
    const ccomplex_t IT_0397 = mtylt_C0iC(9, 0, 0, 0, IT_0003, IT_0005,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0398 = (0 + _Complex_I*2)*M_W*IT_0007*IT_0008;
    const ccomplex_t IT_0399 = -IT_0398;
    const ccomplex_t IT_0400 = IT_0397*IT_0399;
    const ccomplex_t IT_0401 = IT_0020*IT_0177*IT_0400;
    const ccomplex_t IT_0402 = 0.101321183642338*IT_0401;
    const ccomplex_t IT_0403 = mtylt_C0iC(9, 0, 0, 0, IT_0003, IT_0173,
       IT_0017, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0404 = IT_0399*IT_0403;
    const ccomplex_t IT_0405 = IT_0020*IT_0186*IT_0404;
    const ccomplex_t IT_0406 = 0.101321183642338*IT_0405;
    const ccomplex_t IT_0407 = IT_0044*IT_0177*IT_0400;
    const ccomplex_t IT_0408 = 0.101321183642338*IT_0407;
    const ccomplex_t IT_0409 = IT_0044*IT_0186*IT_0404;
    const ccomplex_t IT_0410 = 0.101321183642338*IT_0409;
    const ccomplex_t IT_0411 = 0.101321183642338*m_s;
    const ccomplex_t IT_0412 = mtylt_C0iC(0, 0, 0, 0, IT_0003, IT_0016,
       IT_0017, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0413 = (0 + _Complex_I*2)*IT_0003*IT_0007*IT_0008
      *IT_0032*IT_0034;
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = IT_0020*IT_0067*IT_0412*IT_0414;
    const ccomplex_t IT_0416 = IT_0411*IT_0415;
    const ccomplex_t IT_0417 = 0.101321183642338*m_c;
    const ccomplex_t IT_0418 = mtylt_C0iC(0, 0, 0, 0, IT_0003, IT_0016,
       IT_0005, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0419 = IT_0020*IT_0036*IT_0414*IT_0418;
    const ccomplex_t IT_0420 = IT_0417*IT_0419;
    const ccomplex_t IT_0421 = IT_0044*IT_0074*IT_0414*IT_0418;
    const ccomplex_t IT_0422 = IT_0417*IT_0421;
    const ccomplex_t IT_0423 = IT_0037*IT_0044*IT_0412*IT_0414;
    const ccomplex_t IT_0424 = IT_0411*IT_0423;
    const ccomplex_t IT_0425 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0005, IT_0017,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0426 = (-8)*IT_0425;
    const ccomplex_t IT_0427 = Finite + IT_0426;
    const ccomplex_t IT_0428 = (-4)*IT_0425;
    const ccomplex_t IT_0429 = Finite + IT_0428;
    const ccomplex_t IT_0430 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0005, IT_0017,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0431 = IT_0051*IT_0430;
    const ccomplex_t IT_0432 = 0.5*IT_0427 + -IT_0429 + (-0.5)*IT_0431;
    const ccomplex_t IT_0433 = IT_0425 + IT_0432;
    const ccomplex_t IT_0434 = IT_0010*IT_0082*IT_0121;
    const ccomplex_t IT_0435 = 0.101321183642338*IT_0434;
    const ccomplex_t IT_0436 = IT_0433*IT_0435;
    const ccomplex_t IT_0437 = mtylt_C0iC(0, 0, 0, 0, IT_0005, IT_0173,
       IT_0017, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0438 = IT_0051*IT_0437;
    const ccomplex_t IT_0439 = mtylt_C0iC(9, 0, 0, 0, IT_0005, IT_0173,
       IT_0017, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0440 = (-8)*IT_0439;
    const ccomplex_t IT_0441 = Finite + IT_0440;
    const ccomplex_t IT_0442 = 2*IT_0439 + 0.5*IT_0441;
    const ccomplex_t IT_0443 = IT_0438 + IT_0442;
    const ccomplex_t IT_0444 = IT_0010*IT_0177*IT_0186;
    const ccomplex_t IT_0445 = 0.101321183642338*IT_0444;
    const ccomplex_t IT_0446 = IT_0443*IT_0445;
    const ccomplex_t IT_0447 = IT_0020*IT_0022*IT_0363;
    const ccomplex_t IT_0448 = 0.101321183642338*IT_0447;
    const ccomplex_t IT_0449 = IT_0020*IT_0024*IT_0367;
    const ccomplex_t IT_0450 = 0.101321183642338*IT_0449;
    const ccomplex_t IT_0451 = 0.101321183642338*m_c*m_s;
    const ccomplex_t IT_0452 = mtylt_C0iC(0, 0, 0, 0, IT_0016, IT_0005,
       IT_0017, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0453 = IT_0010*IT_0036*IT_0037*IT_0452;
    const ccomplex_t IT_0454 = IT_0451*IT_0453;
    const ccomplex_t IT_0455 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0003, IT_0005,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0456 = -IT_0081;
    const ccomplex_t IT_0457 = IT_0455*IT_0456;
    const ccomplex_t IT_0458 = IT_0081*IT_0455;
    const ccomplex_t IT_0459 = (-8)*IT_0455;
    const ccomplex_t IT_0460 = Finite + IT_0459;
    const ccomplex_t IT_0461 = 2*IT_0081;
    const ccomplex_t IT_0462 = IT_0460*IT_0461;
    const ccomplex_t IT_0463 = (-2)*IT_0081;
    const ccomplex_t IT_0464 = IT_0460*IT_0463;
    const ccomplex_t IT_0465 = -IT_0458 + (-0.25)*IT_0462 + 0.25*IT_0464;
    const ccomplex_t IT_0466 = IT_0457 + IT_0465;
    const ccomplex_t IT_0467 = IT_0010*IT_0121*IT_0466;
    const ccomplex_t IT_0468 = 0.101321183642338*IT_0467;
    const ccomplex_t IT_0469 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0003, IT_0017,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0470 = IT_0456*IT_0469;
    const ccomplex_t IT_0471 = IT_0081*IT_0469;
    const ccomplex_t IT_0472 = (-8)*IT_0469;
    const ccomplex_t IT_0473 = Finite + IT_0472;
    const ccomplex_t IT_0474 = IT_0461*IT_0473;
    const ccomplex_t IT_0475 = IT_0463*IT_0473;
    const ccomplex_t IT_0476 = -IT_0471 + (-0.25)*IT_0474 + 0.25*IT_0475;
    const ccomplex_t IT_0477 = IT_0470 + IT_0476;
    const ccomplex_t IT_0478 = IT_0010*IT_0082*IT_0477;
    const ccomplex_t IT_0479 = 0.101321183642338*IT_0478;
    const ccomplex_t IT_0480 = (0 + _Complex_I*2)*M_W*IT_0007*IT_0008*IT_0031;
    const ccomplex_t IT_0481 = IT_0360*IT_0480;
    const ccomplex_t IT_0482 = -IT_0480;
    const ccomplex_t IT_0483 = IT_0360*IT_0482;
    const ccomplex_t IT_0484 = (-8)*IT_0360;
    const ccomplex_t IT_0485 = Finite + IT_0484;
    const ccomplex_t IT_0486 = IT_0480*IT_0485;
    const ccomplex_t IT_0487 = IT_0482*IT_0485;
    const ccomplex_t IT_0488 = -IT_0483 + 0.5*IT_0486 + (-0.5)*IT_0487;
    const ccomplex_t IT_0489 = IT_0481 + IT_0488;
    const ccomplex_t IT_0490 = IT_0010*IT_0074*IT_0489;
    const ccomplex_t IT_0491 = 0.101321183642338*IT_0490;
    const ccomplex_t IT_0492 = IT_0366*IT_0480;
    const ccomplex_t IT_0493 = IT_0366*IT_0482;
    const ccomplex_t IT_0494 = (-8)*IT_0366;
    const ccomplex_t IT_0495 = Finite + IT_0494;
    const ccomplex_t IT_0496 = IT_0482*IT_0495;
    const ccomplex_t IT_0497 = IT_0480*IT_0495;
    const ccomplex_t IT_0498 = -IT_0493 + (-0.5)*IT_0496 + 0.5*IT_0497;
    const ccomplex_t IT_0499 = IT_0492 + IT_0498;
    const ccomplex_t IT_0500 = IT_0010*IT_0067*IT_0499;
    const ccomplex_t IT_0501 = 0.101321183642338*IT_0500;
    const ccomplex_t IT_0502 = (-8)*IT_0370;
    const ccomplex_t IT_0503 = Finite + IT_0502;
    const ccomplex_t IT_0504 = IT_0051*IT_0452;
    const ccomplex_t IT_0505 = 0.25*IT_0503 + 0.5*IT_0504;
    const ccomplex_t IT_0506 = IT_0370 + IT_0505;
    const ccomplex_t IT_0507 = IT_0010*IT_0022*IT_0024;
    const ccomplex_t IT_0508 = 0.101321183642338*IT_0507;
    const ccomplex_t IT_0509 = IT_0506*IT_0508;
    const ccomplex_t IT_0510 = IT_0365 + -IT_0369 + 2*IT_0374 + -IT_0379 + 
      -IT_0381 + -IT_0385 + -IT_0387 + -IT_0392 + -IT_0396 + -IT_0402 + IT_0406 
      + -IT_0408 + IT_0410 + -IT_0416 + -IT_0420 + -IT_0422 + -IT_0424 + (-4)
      *IT_0436 + IT_0446 + -IT_0448 + IT_0450 + 2*IT_0454 + (-2)*IT_0468 + 2
      *IT_0479 + 2*IT_0491 + (-2)*IT_0501 + (-2)*IT_0509;
    const ccomplex_t IT_0511 = cpow(2*s_13 + IT_0003 + -IT_0004 + -IT_0005 + 
      -reg_prop, -1);
    const ccomplex_t IT_0512 = IT_0510*IT_0511;
    const ccomplex_t IT_0513 = -IT_0512;
    const ccomplex_t IT_0514 = IT_0012*IT_0047*IT_0513;
    const ccomplex_t IT_0515 = 1.4142135623731*IT_0514;
    const ccomplex_t IT_0516 = (-0.00390625)*IT_0515;
    const ccomplex_t IT_0517 = 0.00390625*IT_0515;
    const ccomplex_t IT_0518 = (-2)*IT_0439 + (-0.5)*IT_0441;
    const ccomplex_t IT_0519 = IT_0438 + IT_0518;
    const ccomplex_t IT_0520 = IT_0445*IT_0519;
    const ccomplex_t IT_0521 = 0.5*IT_0427 + -IT_0429 + 0.5*IT_0431;
    const ccomplex_t IT_0522 = IT_0425 + IT_0521;
    const ccomplex_t IT_0523 = IT_0435*IT_0522;
    const ccomplex_t IT_0524 = 0.25*IT_0503 + (-0.5)*IT_0504;
    const ccomplex_t IT_0525 = IT_0370 + IT_0524;
    const ccomplex_t IT_0526 = IT_0508*IT_0525;
    const ccomplex_t IT_0527 = IT_0365 + IT_0369 + 2*IT_0374 + IT_0379 + 
      -IT_0381 + -IT_0385 + IT_0387 + -IT_0392 + -IT_0396 + IT_0402 + IT_0406 + 
      -IT_0408 + -IT_0410 + -IT_0416 + IT_0420 + -IT_0422 + IT_0424 + IT_0448 +
       IT_0450 + (-2)*IT_0454 + (-2)*IT_0468 + 2*IT_0479 + 2*IT_0491 + (-2)
      *IT_0501 + IT_0520 + (-4)*IT_0523 + 2*IT_0526;
    const ccomplex_t IT_0528 = IT_0511*IT_0527;
    const ccomplex_t IT_0529 = -IT_0528;
    const ccomplex_t IT_0530 = -IT_0529;
    const ccomplex_t IT_0531 = IT_0012*IT_0047*IT_0530;
    const ccomplex_t IT_0532 = 1.4142135623731*IT_0531;
    const ccomplex_t IT_0533 = (-0.00390625)*IT_0532;
    const ccomplex_t IT_0534 = IT_0012*IT_0047*IT_0528;
    const ccomplex_t IT_0535 = 1.4142135623731*IT_0534;
    const ccomplex_t IT_0536 = 0.00390625*IT_0535;
    const ccomplex_t IT_0537 = mtylt_C0iC(9, 0, 0, 0, IT_0003, IT_0004,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0538 = IT_0399*IT_0537;
    const ccomplex_t IT_0539 = IT_0026*IT_0179*IT_0538;
    const ccomplex_t IT_0540 = 0.101321183642338*IT_0539;
    const ccomplex_t IT_0541 = IT_0026*IT_0177*IT_0400;
    const ccomplex_t IT_0542 = 0.101321183642338*IT_0541;
    const ccomplex_t IT_0543 = IT_0097*IT_0179*IT_0538;
    const ccomplex_t IT_0544 = 0.101321183642338*IT_0543;
    const ccomplex_t IT_0545 = mtylt_C0iC(9, 0, 0, 0, IT_0003, IT_0016,
       IT_0004, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0546 = IT_0480*IT_0545;
    const ccomplex_t IT_0547 = (-8)*IT_0545;
    const ccomplex_t IT_0548 = Finite + IT_0547;
    const ccomplex_t IT_0549 = IT_0480*IT_0548;
    const ccomplex_t IT_0550 = IT_0482*IT_0545;
    const ccomplex_t IT_0551 = IT_0482*IT_0548;
    const ccomplex_t IT_0552 = 0.5*IT_0549 + -IT_0550 + (-0.5)*IT_0551;
    const ccomplex_t IT_0553 = IT_0546 + IT_0552;
    const ccomplex_t IT_0554 = IT_0012*IT_0067*IT_0553;
    const ccomplex_t IT_0555 = 0.101321183642338*IT_0554;
    const ccomplex_t IT_0556 = IT_0012*IT_0177*IT_0389*IT_0390;
    const ccomplex_t IT_0557 = 0.101321183642338*IT_0556;
    const ccomplex_t IT_0558 = IT_0455*IT_0461;
    const ccomplex_t IT_0559 = IT_0456*IT_0460;
    const ccomplex_t IT_0560 = 0.25*IT_0462 + (-0.5)*IT_0559;
    const ccomplex_t IT_0561 = IT_0558 + IT_0560;
    const ccomplex_t IT_0562 = IT_0012*IT_0121*IT_0561;
    const ccomplex_t IT_0563 = 0.101321183642338*IT_0562;
    const ccomplex_t IT_0564 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0004, IT_0005,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0565 = (-8)*IT_0564;
    const ccomplex_t IT_0566 = Finite + IT_0565;
    const ccomplex_t IT_0567 = (-4)*IT_0564;
    const ccomplex_t IT_0568 = Finite + IT_0567;
    const ccomplex_t IT_0569 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0004, IT_0005,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0570 = IT_0136*IT_0569;
    const ccomplex_t IT_0571 = 0.5*IT_0566 + -IT_0568 + (-0.5)*IT_0570;
    const ccomplex_t IT_0572 = IT_0564 + IT_0571;
    const ccomplex_t IT_0573 = IT_0012*IT_0082*IT_0121;
    const ccomplex_t IT_0574 = 0.101321183642338*IT_0573;
    const ccomplex_t IT_0575 = IT_0572*IT_0574;
    const ccomplex_t IT_0576 = IT_0097*IT_0121*IT_0376*IT_0377;
    const ccomplex_t IT_0577 = 0.101321183642338*IT_0576;
    const ccomplex_t IT_0578 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0003, IT_0004,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0579 = m_b*IT_0578;
    const ccomplex_t IT_0580 = IT_0026*IT_0082*IT_0377*IT_0579;
    const ccomplex_t IT_0581 = 0.101321183642338*IT_0580;
    const ccomplex_t IT_0582 = mtylt_C0iC(0, 0, 0, 0, IT_0003, IT_0004,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0583 = m_b*IT_0582;
    const ccomplex_t IT_0584 = IT_0012*IT_0179*IT_0390*IT_0583;
    const ccomplex_t IT_0585 = 0.101321183642338*IT_0584;
    const ccomplex_t IT_0586 = IT_0036*IT_0097*IT_0414*IT_0418;
    const ccomplex_t IT_0587 = IT_0417*IT_0586;
    const ccomplex_t IT_0588 = 0.101321183642338*m_b*m_c;
    const ccomplex_t IT_0589 = mtylt_C0iC(0, 0, 0, 0, IT_0016, IT_0004,
       IT_0005, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0590 = IT_0012*IT_0036*IT_0037*IT_0589;
    const ccomplex_t IT_0591 = IT_0588*IT_0590;
    const ccomplex_t IT_0592 = IT_0012*IT_0074*IT_0489;
    const ccomplex_t IT_0593 = 0.101321183642338*IT_0592;
    const ccomplex_t IT_0594 = IT_0361*IT_0545;
    const ccomplex_t IT_0595 = IT_0097*IT_0192*IT_0594;
    const ccomplex_t IT_0596 = 0.101321183642338*IT_0595;
    const ccomplex_t IT_0597 = IT_0026*IT_0121*IT_0376*IT_0377;
    const ccomplex_t IT_0598 = 0.101321183642338*IT_0597;
    const ccomplex_t IT_0599 = IT_0082*IT_0097*IT_0377*IT_0579;
    const ccomplex_t IT_0600 = 0.101321183642338*IT_0599;
    const ccomplex_t IT_0601 = IT_0097*IT_0177*IT_0400;
    const ccomplex_t IT_0602 = 0.101321183642338*IT_0601;
    const ccomplex_t IT_0603 = IT_0026*IT_0074*IT_0414*IT_0418;
    const ccomplex_t IT_0604 = IT_0417*IT_0603;
    const ccomplex_t IT_0605 = 0.101321183642338*m_b;
    const ccomplex_t IT_0606 = mtylt_C0iC(0, 0, 0, 0, IT_0003, IT_0016,
       IT_0004, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0607 = IT_0026*IT_0037*IT_0414*IT_0606;
    const ccomplex_t IT_0608 = IT_0605*IT_0607;
    const ccomplex_t IT_0609 = IT_0067*IT_0097*IT_0414*IT_0606;
    const ccomplex_t IT_0610 = IT_0605*IT_0609;
    const ccomplex_t IT_0611 = IT_0026*IT_0192*IT_0594;
    const ccomplex_t IT_0612 = 0.101321183642338*IT_0611;
    const ccomplex_t IT_0613 = IT_0360*IT_0361;
    const ccomplex_t IT_0614 = IT_0022*IT_0026*IT_0613;
    const ccomplex_t IT_0615 = 0.101321183642338*IT_0614;
    const ccomplex_t IT_0616 = IT_0022*IT_0097*IT_0613;
    const ccomplex_t IT_0617 = 0.101321183642338*IT_0616;
    const ccomplex_t IT_0618 = mtylt_C0iC(9, 0, 0, 0, IT_0016, IT_0004,
       IT_0005, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0619 = (-2)*IT_0618;
    const ccomplex_t IT_0620 = Finite + IT_0619;
    const ccomplex_t IT_0621 = IT_0012*IT_0067*IT_0074*IT_0620;
    const ccomplex_t IT_0622 = 0.101321183642338*IT_0621;
    const ccomplex_t IT_0623 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0003, IT_0004,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0624 = IT_0461*IT_0623;
    const ccomplex_t IT_0625 = (-8)*IT_0623;
    const ccomplex_t IT_0626 = Finite + IT_0625;
    const ccomplex_t IT_0627 = IT_0456*IT_0626;
    const ccomplex_t IT_0628 = IT_0461*IT_0626;
    const ccomplex_t IT_0629 = (-0.5)*IT_0627 + 0.25*IT_0628;
    const ccomplex_t IT_0630 = IT_0624 + IT_0629;
    const ccomplex_t IT_0631 = IT_0012*IT_0082*IT_0630;
    const ccomplex_t IT_0632 = 0.101321183642338*IT_0631;
    const ccomplex_t IT_0633 = mtylt_C0iC(0, 0, 0, 0, IT_0004, IT_0005,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0634 = IT_0136*IT_0633;
    const ccomplex_t IT_0635 = mtylt_C0iC(9, 0, 0, 0, IT_0004, IT_0005,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0636 = (-8)*IT_0635;
    const ccomplex_t IT_0637 = Finite + IT_0636;
    const ccomplex_t IT_0638 = 2*IT_0635 + 0.5*IT_0637;
    const ccomplex_t IT_0639 = IT_0634 + IT_0638;
    const ccomplex_t IT_0640 = IT_0012*IT_0177*IT_0179;
    const ccomplex_t IT_0641 = 0.101321183642338*IT_0640;
    const ccomplex_t IT_0642 = IT_0639*IT_0641;
    const ccomplex_t IT_0643 = IT_0136*IT_0589;
    const ccomplex_t IT_0644 = (-8)*IT_0618;
    const ccomplex_t IT_0645 = Finite + IT_0644;
    const ccomplex_t IT_0646 = 0.5*IT_0643 + 0.25*IT_0645;
    const ccomplex_t IT_0647 = IT_0618 + IT_0646;
    const ccomplex_t IT_0648 = IT_0012*IT_0022*IT_0192;
    const ccomplex_t IT_0649 = 0.101321183642338*IT_0648;
    const ccomplex_t IT_0650 = IT_0647*IT_0649;
    const ccomplex_t IT_0651 = IT_0540 + -IT_0542 + IT_0544 + (-2)*IT_0555 + 
      -IT_0557 + 2*IT_0563 + (-4)*IT_0575 + -IT_0577 + -IT_0581 + -IT_0585 + 
      -IT_0587 + 2*IT_0591 + 2*IT_0593 + -IT_0596 + -IT_0598 + -IT_0600 + 
      -IT_0602 + -IT_0604 + -IT_0608 + -IT_0610 + IT_0612 + -IT_0615 + IT_0617 +
       2*IT_0622 + (-2)*IT_0632 + IT_0642 + (-2)*IT_0650;
    const ccomplex_t IT_0652 = cpow((-2)*s_24 + IT_0003 + -IT_0005 + -IT_0017 
      + -reg_prop, -1);
    const ccomplex_t IT_0653 = IT_0651*IT_0652;
    const ccomplex_t IT_0654 = -IT_0653;
    const ccomplex_t IT_0655 = IT_0010*IT_0047*IT_0654;
    const ccomplex_t IT_0656 = 1.4142135623731*IT_0655;
    const ccomplex_t IT_0657 = (-0.00390625)*IT_0656;
    const ccomplex_t IT_0658 = 0.5*IT_0566 + -IT_0568 + 0.5*IT_0570;
    const ccomplex_t IT_0659 = IT_0564 + IT_0658;
    const ccomplex_t IT_0660 = IT_0574*IT_0659;
    const ccomplex_t IT_0661 = (-2)*IT_0635 + (-0.5)*IT_0637;
    const ccomplex_t IT_0662 = IT_0634 + IT_0661;
    const ccomplex_t IT_0663 = IT_0641*IT_0662;
    const ccomplex_t IT_0664 = (-0.5)*IT_0643 + 0.25*IT_0645;
    const ccomplex_t IT_0665 = IT_0618 + IT_0664;
    const ccomplex_t IT_0666 = IT_0649*IT_0665;
    const ccomplex_t IT_0667 = IT_0540 + IT_0542 + -IT_0544 + 2*IT_0555 +
       IT_0557 + (-2)*IT_0563 + -IT_0577 + -IT_0581 + IT_0585 + -IT_0587 + 2
      *IT_0591 + (-2)*IT_0593 + IT_0596 + IT_0598 + IT_0600 + -IT_0602 + IT_0604
       + -IT_0608 + IT_0610 + IT_0612 + IT_0615 + IT_0617 + (-2)*IT_0622 + 2
      *IT_0632 + 4*IT_0660 + -IT_0663 + (-2)*IT_0666;
    const ccomplex_t IT_0668 = IT_0652*IT_0667;
    const ccomplex_t IT_0669 = IT_0010*IT_0047*IT_0668;
    const ccomplex_t IT_0670 = 1.4142135623731*IT_0669;
    const ccomplex_t IT_0671 = 0.00390625*IT_0670;
    const ccomplex_t IT_0672 = 0.00390625*IT_0656;
    const ccomplex_t IT_0673 = (-0.00390625)*IT_0670;
    const ccomplex_t IT_0674 = (0 + _Complex_I*0.353553390593274)*IT_0000
      *IT_0001*IT_0002*IT_0015 + (0 + _Complex_I*(-0.0055242717280199))*IT_0000
      *IT_0001*IT_0002*IT_0018*IT_0028 + (0 + _Complex_I*(-0.0055242717280199))
      *IT_0000*IT_0001*IT_0002*IT_0039 + (0 + _Complex_I*(-0.0055242717280199))
      *IT_0000*IT_0001*IT_0002*IT_0046 + (0 + _Complex_I*0.00390625)*IT_0050 + 
      (0 + _Complex_I*(-0.0055242717280199))*IT_0000*IT_0001*IT_0002*IT_0053
      *IT_0055 + (0 + _Complex_I*0.0055242717280199)*IT_0000*IT_0001*IT_0002
      *IT_0018*IT_0057 + (0 + _Complex_I*1)*IT_0220 + (0 + _Complex_I*1)*IT_0323
       + (0 + _Complex_I*-1)*IT_0329 + (0 + _Complex_I*-1)*IT_0337 + (0 +
       _Complex_I*-1)*IT_0342 + (0 + _Complex_I*-1)*IT_0348 + (0 + _Complex_I*1)
      *IT_0353 + (0 + _Complex_I*1)*IT_0359 + IT_0516 + -IT_0517 + -IT_0533 +
       IT_0536 + IT_0657 + -IT_0671 + -IT_0672 + IT_0673;
    return create_ccomplex_return(IT_0674);
}

