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
    const ccomplex_t IT_0003 = (0 + _Complex_I*0.101321183642338)*m_b*m_c;
    const ccomplex_t IT_0004 = pow(M_W, 2);
    const ccomplex_t IT_0005 = pow(M_Z, 2);
    const ccomplex_t IT_0006 = pow(m_b, 2);
    const ccomplex_t IT_0007 = pow(m_c, 2);
    const ccomplex_t IT_0008 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0005, IT_0006, IT_0007);
    const ccomplex_t IT_0009 = pow(1.4142135623731, 0.5);
    const ccomplex_t IT_0010 = pow(G_F, 0.5);
    const ccomplex_t IT_0011 = (0 + _Complex_I*2.82842712474619)*m_s*conj(V_cs
      )*IT_0009*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*2.82842712474619)*m_c*V_cb
      *IT_0009*IT_0010;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, 2);
    const ccomplex_t IT_0019 = (0 + _Complex_I*2)*M_W*IT_0009*IT_0010*IT_0016
      *IT_0018;
    const ccomplex_t IT_0020 = (-0.666666666666667)*IT_0019;
    const ccomplex_t IT_0021 = 0.333333333333333*IT_0019;
    const ccomplex_t IT_0026 = pow(m_s, 2);
    const ccomplex_t IT_0027 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0005, IT_0006, IT_0026);
    const ccomplex_t IT_0028 = (0 + _Complex_I*2.82842712474619)*m_c*conj(V_cs
      )*IT_0009*IT_0010;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0032 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0005, IT_0007, IT_0026);
    const ccomplex_t IT_0033 = 2*m_c*IT_0009*IT_0010;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = 2*m_s*IT_0009*IT_0010;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0014*IT_0029*IT_0034*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*0.101321183642338)*IT_0037;
    const ccomplex_t IT_0039 = IT_0012*IT_0014*IT_0034*IT_0036;
    const ccomplex_t IT_0040 = (0 + _Complex_I*0.101321183642338)*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*2.82842712474619)*m_b*V_cb
      *IT_0009*IT_0010;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0012*IT_0034*IT_0036*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*0.101321183642338)*IT_0043;
    const ccomplex_t IT_0045 = IT_0000*IT_0001*IT_0002;
    const ccomplex_t IT_0046 = IT_0016*IT_0017;
    const ccomplex_t IT_0047 = IT_0009*IT_0010*IT_0017;
    const ccomplex_t IT_0048 = 4*M_W*IT_0046*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = cpow(IT_0017, -1);
    const ccomplex_t IT_0051 = IT_0015*IT_0050;
    const ccomplex_t IT_0052 = 4*M_W*IT_0047*IT_0051;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0049 + 3*IT_0053);
    const ccomplex_t IT_0055 = (-0.166666666666667)*IT_0054;
    const ccomplex_t IT_0056 = cpow(IT_0055, 2);
    const ccomplex_t IT_0057 = (0 + _Complex_I*0.101321183642338)*IT_0056;
    const ccomplex_t IT_0058 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0005, IT_0006, IT_0026);
    const ccomplex_t IT_0059 = (0 + _Complex_I*2.82842712474619)*M_W*conj(V_cs
      )*IT_0009*IT_0010;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*2.82842712474619)*M_W*V_cb
      *IT_0009*IT_0010;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = IT_0058*IT_0060*IT_0062;
    const ccomplex_t IT_0064 = IT_0057*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0049 + (-3)*IT_0053);
    const ccomplex_t IT_0066 = (-0.166666666666667)*IT_0065;
    const ccomplex_t IT_0067 = cpow(IT_0066, 2);
    const ccomplex_t IT_0068 = (0 + _Complex_I*0.101321183642338)*IT_0067;
    const ccomplex_t IT_0069 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0005, IT_0007, IT_0007);
    const ccomplex_t IT_0070 = IT_0060*IT_0062*IT_0069;
    const ccomplex_t IT_0071 = IT_0068*IT_0070;
    const ccomplex_t IT_0072 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0006, IT_0026);
    const ccomplex_t IT_0073 = IT_0060*IT_0062;
    const ccomplex_t IT_0074 = (0 + _Complex_I*2)*M_W*IT_0009*IT_0010*IT_0017;
    const ccomplex_t IT_0075 = (-0.333333333333333)*IT_0074;
    const ccomplex_t IT_0076 = cpow(IT_0075, 2);
    const ccomplex_t IT_0077 = (0 + _Complex_I*0.101321183642338)*IT_0076;
    const ccomplex_t IT_0078 = IT_0073*IT_0077;
    const ccomplex_t IT_0079 = IT_0072*IT_0078;
    const ccomplex_t IT_0080 = IT_0064 + IT_0071 + IT_0079;
    const ccomplex_t IT_0081 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0005, IT_0007, IT_0007);
    const ccomplex_t IT_0082 = IT_0014*IT_0029*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*0.101321183642338)*IT_0007
      *IT_0067;
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = IT_0032*IT_0055*IT_0060*IT_0062*IT_0066;
    const ccomplex_t IT_0086 = (0 + _Complex_I*0.101321183642338)*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*0.101321183642338)*m_c*m_s;
    const ccomplex_t IT_0088 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0005, IT_0007, IT_0026);
    const ccomplex_t IT_0089 = IT_0012*IT_0014*IT_0055*IT_0066*IT_0088;
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0005, IT_0006, IT_0007);
    const ccomplex_t IT_0092 = IT_0055*IT_0060*IT_0062*IT_0066*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*0.101321183642338)*IT_0092;
    const ccomplex_t IT_0094 = IT_0008*IT_0029*IT_0042*IT_0055*IT_0066;
    const ccomplex_t IT_0095 = IT_0003*IT_0094;
    const ccomplex_t IT_0096 = IT_0008*IT_0014*IT_0021*IT_0029*IT_0066;
    const ccomplex_t IT_0097 = IT_0003*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*0.101321183642338)*IT_0007;
    const ccomplex_t IT_0099 = IT_0012*IT_0014*IT_0020*IT_0066*IT_0081;
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = IT_0014*IT_0021*IT_0029*IT_0066*IT_0088;
    const ccomplex_t IT_0102 = IT_0087*IT_0101;
    const ccomplex_t IT_0103 = IT_0020*IT_0029*IT_0042*IT_0066*IT_0081;
    const ccomplex_t IT_0104 = IT_0098*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*0.101321183642338)*m_b*m_s;
    const ccomplex_t IT_0106 = IT_0012*IT_0014*IT_0021*IT_0027*IT_0055;
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = IT_0012*IT_0020*IT_0042*IT_0055*IT_0088;
    const ccomplex_t IT_0109 = IT_0087*IT_0108;
    const ccomplex_t IT_0110 = IT_0008*IT_0012*IT_0020*IT_0042*IT_0055;
    const ccomplex_t IT_0111 = IT_0003*IT_0110;
    const ccomplex_t IT_0112 = IT_0021*IT_0027*IT_0029*IT_0042*IT_0055;
    const ccomplex_t IT_0113 = IT_0105*IT_0112;
    const ccomplex_t IT_0114 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0007, IT_0007);
    const ccomplex_t IT_0115 = IT_0007*IT_0114;
    const ccomplex_t IT_0116 = 0.666666666666667*IT_0074;
    const ccomplex_t IT_0117 = cpow(IT_0116, 2);
    const ccomplex_t IT_0118 = (0 + _Complex_I*0.101321183642338)*IT_0117;
    const ccomplex_t IT_0119 = IT_0029*IT_0042;
    const ccomplex_t IT_0120 = IT_0118*IT_0119;
    const ccomplex_t IT_0121 = IT_0115*IT_0120;
    const ccomplex_t IT_0122 = IT_0012*IT_0042;
    const ccomplex_t IT_0123 = IT_0118*IT_0122;
    const ccomplex_t IT_0124 = IT_0115*IT_0123;
    const ccomplex_t IT_0125 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0007, IT_0026);
    const ccomplex_t IT_0126 = m_c*m_s;
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = IT_0014*IT_0029*IT_0075*IT_0116;
    const ccomplex_t IT_0129 = (0 + _Complex_I*0.101321183642338)*IT_0128;
    const ccomplex_t IT_0130 = IT_0127*IT_0129;
    const ccomplex_t IT_0131 = IT_0029*IT_0042*IT_0075*IT_0116;
    const ccomplex_t IT_0132 = (0 + _Complex_I*0.101321183642338)*IT_0131;
    const ccomplex_t IT_0133 = IT_0127*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*0.101321183642338)*m_b*m_s
      *IT_0056;
    const ccomplex_t IT_0135 = IT_0012*IT_0027*IT_0042;
    const ccomplex_t IT_0136 = IT_0134*IT_0135;
    const ccomplex_t IT_0137 = IT_0012*IT_0014;
    const ccomplex_t IT_0138 = IT_0118*IT_0137;
    const ccomplex_t IT_0139 = IT_0115*IT_0138;
    const ccomplex_t IT_0140 = IT_0014*IT_0029;
    const ccomplex_t IT_0141 = IT_0118*IT_0140;
    const ccomplex_t IT_0142 = IT_0115*IT_0141;
    const ccomplex_t IT_0143 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0007, IT_0026);
    const ccomplex_t IT_0144 = IT_0060*IT_0062*IT_0075*IT_0116;
    const ccomplex_t IT_0145 = (0 + _Complex_I*0.101321183642338)*IT_0144;
    const ccomplex_t IT_0146 = IT_0143*IT_0145;
    const ccomplex_t IT_0147 = IT_0012*IT_0014*IT_0075*IT_0116;
    const ccomplex_t IT_0148 = (0 + _Complex_I*0.101321183642338)*IT_0147;
    const ccomplex_t IT_0149 = IT_0127*IT_0148;
    const ccomplex_t IT_0150 = IT_0012*IT_0042*IT_0075*IT_0116;
    const ccomplex_t IT_0151 = (0 + _Complex_I*0.101321183642338)*IT_0150;
    const ccomplex_t IT_0152 = IT_0127*IT_0151;
    const ccomplex_t IT_0153 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0006, IT_0007);
    const ccomplex_t IT_0154 = m_b*m_c;
    const ccomplex_t IT_0155 = IT_0153*IT_0154;
    const ccomplex_t IT_0156 = IT_0148*IT_0155;
    const ccomplex_t IT_0157 = IT_0129*IT_0155;
    const ccomplex_t IT_0158 = IT_0151*IT_0155;
    const ccomplex_t IT_0159 = IT_0132*IT_0155;
    const ccomplex_t IT_0160 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0006, IT_0007);
    const ccomplex_t IT_0161 = IT_0145*IT_0160;
    const ccomplex_t IT_0162 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0006, IT_0026);
    const ccomplex_t IT_0163 = m_b*m_s;
    const ccomplex_t IT_0164 = IT_0162*IT_0163;
    const ccomplex_t IT_0165 = IT_0077*IT_0137;
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = IT_0077*IT_0140;
    const ccomplex_t IT_0168 = IT_0164*IT_0167;
    const ccomplex_t IT_0169 = IT_0077*IT_0122;
    const ccomplex_t IT_0170 = IT_0164*IT_0169;
    const ccomplex_t IT_0171 = IT_0077*IT_0119;
    const ccomplex_t IT_0172 = IT_0164*IT_0171;
    const ccomplex_t IT_0173 = pow(m_h, 2);
    const ccomplex_t IT_0174 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0006, IT_0007, IT_0173);
    const ccomplex_t IT_0175 = IT_0154*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*2)*m_c*IT_0009*IT_0010;
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*2)*m_b*IT_0009*IT_0010;
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = IT_0060*IT_0062*IT_0177*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*0.101321183642338)*IT_0180;
    const ccomplex_t IT_0182 = IT_0175*IT_0181;
    const ccomplex_t IT_0183 = IT_0008*IT_0154;
    const ccomplex_t IT_0184 = 2*m_b*IT_0009*IT_0010;
    const ccomplex_t IT_0185 = 0.5*IT_0184;
    const ccomplex_t IT_0186 = IT_0034*IT_0060*IT_0062*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*0.101321183642338)*IT_0186;
    const ccomplex_t IT_0188 = IT_0183*IT_0187;
    const ccomplex_t IT_0189 = IT_0027*IT_0163;
    const ccomplex_t IT_0190 = IT_0036*IT_0060*IT_0062*IT_0185;
    const ccomplex_t IT_0191 = (0 + _Complex_I*0.101321183642338)*IT_0190;
    const ccomplex_t IT_0192 = IT_0189*IT_0191;
    const ccomplex_t IT_0193 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0007, IT_0007, IT_0173);
    const ccomplex_t IT_0194 = IT_0007*IT_0193;
    const ccomplex_t IT_0195 = cpow(IT_0177, 2);
    const ccomplex_t IT_0196 = (0 + _Complex_I*0.101321183642338)*IT_0195;
    const ccomplex_t IT_0197 = IT_0073*IT_0196;
    const ccomplex_t IT_0198 = IT_0194*IT_0197;
    const ccomplex_t IT_0199 = IT_0007*IT_0081;
    const ccomplex_t IT_0200 = cpow(IT_0034, 2);
    const ccomplex_t IT_0201 = (0 + _Complex_I*0.101321183642338)*IT_0200;
    const ccomplex_t IT_0202 = IT_0073*IT_0201;
    const ccomplex_t IT_0203 = IT_0199*IT_0202;
    const ccomplex_t IT_0204 = IT_0088*IT_0126;
    const ccomplex_t IT_0205 = IT_0034*IT_0036*IT_0060*IT_0062;
    const ccomplex_t IT_0206 = (0 + _Complex_I*0.101321183642338)*IT_0205;
    const ccomplex_t IT_0207 = IT_0204*IT_0206;
    const ccomplex_t IT_0208 = IT_0032*IT_0044;
    const ccomplex_t IT_0209 = IT_0029*IT_0034*IT_0036*IT_0042;
    const ccomplex_t IT_0210 = (0 + _Complex_I*0.101321183642338)*IT_0209;
    const ccomplex_t IT_0211 = IT_0032*IT_0210;
    const ccomplex_t IT_0212 = (-0.25)*IT_0084 + (-4)*IT_0086 + (-0.25)
      *IT_0090 + (-4)*IT_0093 + (-0.25)*IT_0095 + (-0.25)*IT_0097 + (-0.25)
      *IT_0100 + (-0.25)*IT_0102 + (-0.25)*IT_0104 + (-0.25)*IT_0107 + (-0.25)
      *IT_0109 + (-0.25)*IT_0111 + (-0.25)*IT_0113 + (-0.25)*IT_0121 + (-0.25)
      *IT_0124 + (-0.25)*IT_0130 + (-0.25)*IT_0133 + (-0.25)*IT_0136 + (-0.25)
      *IT_0139 + (-0.25)*IT_0142 + (-4)*IT_0146 + (-0.25)*IT_0149 + (-0.25)
      *IT_0152 + (-0.25)*IT_0156 + (-0.25)*IT_0157 + (-0.25)*IT_0158 + (-0.25)
      *IT_0159 + (-4)*IT_0161 + (-0.25)*IT_0166 + (-0.25)*IT_0168 + (-0.25)
      *IT_0170 + (-0.25)*IT_0172 + (-0.25)*IT_0182 + (-0.25)*IT_0188 + 0.25
      *IT_0192 + (-0.25)*IT_0198 + 0.25*IT_0203 + (-0.25)*IT_0207 + 0.25*IT_0208
       + (-0.25)*IT_0211;
    const ccomplex_t IT_0213 = IT_0080 + IT_0212;
    const ccomplex_t IT_0214 = IT_0045*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = 0.015625*IT_0215;
    const ccomplex_t IT_0217 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0007, IT_0007);
    const ccomplex_t IT_0218 = IT_0073*IT_0118;
    const ccomplex_t IT_0219 = IT_0217*IT_0218;
    const ccomplex_t IT_0220 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0006, IT_0007, IT_0173);
    const ccomplex_t IT_0221 = IT_0012*IT_0014*IT_0177*IT_0179;
    const ccomplex_t IT_0222 = (0 + _Complex_I*0.101321183642338)*IT_0221;
    const ccomplex_t IT_0223 = IT_0220*IT_0222;
    const ccomplex_t IT_0224 = IT_0012*IT_0034*IT_0042*IT_0185;
    const ccomplex_t IT_0225 = (0 + _Complex_I*0.101321183642338)*IT_0224;
    const ccomplex_t IT_0226 = IT_0091*IT_0225;
    const ccomplex_t IT_0227 = IT_0012*IT_0036*IT_0042*IT_0185;
    const ccomplex_t IT_0228 = (0 + _Complex_I*0.101321183642338)*IT_0227;
    const ccomplex_t IT_0229 = IT_0058*IT_0228;
    const ccomplex_t IT_0230 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0007, IT_0007, IT_0173);
    const ccomplex_t IT_0231 = IT_0140*IT_0196;
    const ccomplex_t IT_0232 = IT_0230*IT_0231;
    const ccomplex_t IT_0233 = IT_0122*IT_0196;
    const ccomplex_t IT_0234 = IT_0230*IT_0233;
    const ccomplex_t IT_0235 = IT_0119*IT_0196;
    const ccomplex_t IT_0236 = IT_0230*IT_0235;
    const ccomplex_t IT_0237 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0007, IT_0173, IT_0026);
    const ccomplex_t IT_0238 = (0 + _Complex_I*2)*m_s*IT_0009*IT_0010;
    const ccomplex_t IT_0239 = (-0.5)*IT_0238;
    const ccomplex_t IT_0240 = IT_0029*IT_0042*IT_0177*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*0.101321183642338)*IT_0240;
    const ccomplex_t IT_0242 = IT_0237*IT_0241;
    const ccomplex_t IT_0243 = IT_0029*IT_0034*IT_0042*IT_0185;
    const ccomplex_t IT_0244 = (0 + _Complex_I*0.101321183642338)*IT_0243;
    const ccomplex_t IT_0245 = IT_0091*IT_0244;
    const ccomplex_t IT_0246 = IT_0012*IT_0014*IT_0036*IT_0185;
    const ccomplex_t IT_0247 = (0 + _Complex_I*0.101321183642338)*IT_0246;
    const ccomplex_t IT_0248 = IT_0058*IT_0247;
    const ccomplex_t IT_0249 = IT_0122*IT_0201;
    const ccomplex_t IT_0250 = IT_0069*IT_0249;
    const ccomplex_t IT_0251 = mtylt_D0iC(12, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0006, IT_0173, IT_0026);
    const ccomplex_t IT_0252 = IT_0012*IT_0042*IT_0179*IT_0239;
    const ccomplex_t IT_0253 = (0 + _Complex_I*0.101321183642338)*IT_0252;
    const ccomplex_t IT_0254 = IT_0251*IT_0253;
    const ccomplex_t IT_0255 = IT_0014*IT_0029*IT_0036*IT_0185;
    const ccomplex_t IT_0256 = (0 + _Complex_I*0.101321183642338)*IT_0255;
    const ccomplex_t IT_0257 = IT_0058*IT_0256;
    const ccomplex_t IT_0258 = IT_0137*IT_0196;
    const ccomplex_t IT_0259 = IT_0230*IT_0258;
    const ccomplex_t IT_0260 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0007, IT_0173, IT_0026);
    const ccomplex_t IT_0261 = IT_0126*IT_0260;
    const ccomplex_t IT_0262 = IT_0060*IT_0062*IT_0177*IT_0239;
    const ccomplex_t IT_0263 = (0 + _Complex_I*0.101321183642338)*IT_0262;
    const ccomplex_t IT_0264 = IT_0261*IT_0263;
    const ccomplex_t IT_0265 = IT_0014*IT_0029*IT_0177*IT_0239;
    const ccomplex_t IT_0266 = (0 + _Complex_I*0.101321183642338)*IT_0265;
    const ccomplex_t IT_0267 = IT_0237*IT_0266;
    const ccomplex_t IT_0268 = IT_0014*IT_0029*IT_0177*IT_0179;
    const ccomplex_t IT_0269 = (0 + _Complex_I*0.101321183642338)*IT_0268;
    const ccomplex_t IT_0270 = IT_0220*IT_0269;
    const ccomplex_t IT_0271 = IT_0012*IT_0042*IT_0177*IT_0179;
    const ccomplex_t IT_0272 = (0 + _Complex_I*0.101321183642338)*IT_0271;
    const ccomplex_t IT_0273 = IT_0220*IT_0272;
    const ccomplex_t IT_0274 = IT_0029*IT_0042*IT_0177*IT_0179;
    const ccomplex_t IT_0275 = (0 + _Complex_I*0.101321183642338)*IT_0274;
    const ccomplex_t IT_0276 = IT_0220*IT_0275;
    const ccomplex_t IT_0277 = mtylt_D0iC(0, 0, 0, 0, 0, 0, 0, IT_0004,
       IT_0006, IT_0173, IT_0026);
    const ccomplex_t IT_0278 = IT_0163*IT_0277;
    const ccomplex_t IT_0279 = IT_0060*IT_0062*IT_0179*IT_0239;
    const ccomplex_t IT_0280 = (0 + _Complex_I*0.101321183642338)*IT_0279;
    const ccomplex_t IT_0281 = IT_0278*IT_0280;
    const ccomplex_t IT_0282 = IT_0012*IT_0014*IT_0179*IT_0239;
    const ccomplex_t IT_0283 = (0 + _Complex_I*0.101321183642338)*IT_0282;
    const ccomplex_t IT_0284 = IT_0251*IT_0283;
    const ccomplex_t IT_0285 = IT_0014*IT_0029*IT_0179*IT_0239;
    const ccomplex_t IT_0286 = (0 + _Complex_I*0.101321183642338)*IT_0285;
    const ccomplex_t IT_0287 = IT_0251*IT_0286;
    const ccomplex_t IT_0288 = IT_0029*IT_0042*IT_0179*IT_0239;
    const ccomplex_t IT_0289 = (0 + _Complex_I*0.101321183642338)*IT_0288;
    const ccomplex_t IT_0290 = IT_0251*IT_0289;
    const ccomplex_t IT_0291 = IT_0012*IT_0014*IT_0034*IT_0185;
    const ccomplex_t IT_0292 = (0 + _Complex_I*0.101321183642338)*IT_0291;
    const ccomplex_t IT_0293 = IT_0091*IT_0292;
    const ccomplex_t IT_0294 = IT_0014*IT_0029*IT_0034*IT_0185;
    const ccomplex_t IT_0295 = (0 + _Complex_I*0.101321183642338)*IT_0294;
    const ccomplex_t IT_0296 = IT_0091*IT_0295;
    const ccomplex_t IT_0297 = IT_0029*IT_0036*IT_0042*IT_0185;
    const ccomplex_t IT_0298 = (0 + _Complex_I*0.101321183642338)*IT_0297;
    const ccomplex_t IT_0299 = IT_0058*IT_0298;
    const ccomplex_t IT_0300 = IT_0012*IT_0014*IT_0177*IT_0239;
    const ccomplex_t IT_0301 = (0 + _Complex_I*0.101321183642338)*IT_0300;
    const ccomplex_t IT_0302 = IT_0237*IT_0301;
    const ccomplex_t IT_0303 = IT_0012*IT_0042*IT_0177*IT_0239;
    const ccomplex_t IT_0304 = (0 + _Complex_I*0.101321183642338)*IT_0303;
    const ccomplex_t IT_0305 = IT_0237*IT_0304;
    const ccomplex_t IT_0306 = IT_0137*IT_0201;
    const ccomplex_t IT_0307 = IT_0069*IT_0306;
    const ccomplex_t IT_0308 = IT_0140*IT_0201;
    const ccomplex_t IT_0309 = IT_0069*IT_0308;
    const ccomplex_t IT_0310 = IT_0119*IT_0201;
    const ccomplex_t IT_0311 = IT_0069*IT_0310;
    const ccomplex_t IT_0312 = (-0.25)*IT_0223 + 0.25*IT_0226 + (-0.25)
      *IT_0229 + 0.25*IT_0232 + 0.25*IT_0234 + 0.25*IT_0236 + (-0.25)*IT_0242 + 
      (-0.25)*IT_0245 + 0.25*IT_0248 + (-0.25)*IT_0250 + 0.25*IT_0254 + (-0.25)
      *IT_0257 + 0.25*IT_0259 + (-0.25)*IT_0264 + (-0.25)*IT_0267 + (-0.25)
      *IT_0270 + (-0.25)*IT_0273 + (-0.25)*IT_0276 + (-0.25)*IT_0281 + 0.25
      *IT_0284 + 0.25*IT_0287 + 0.25*IT_0290 + (-0.25)*IT_0293 + 0.25*IT_0296 +
       0.25*IT_0299 + (-0.25)*IT_0302 + (-0.25)*IT_0305 + 0.25*IT_0307 + (-0.25)
      *IT_0309 + 0.25*IT_0311;
    const ccomplex_t IT_0313 = IT_0219 + IT_0312;
    const ccomplex_t IT_0314 = IT_0045*IT_0313;
    const ccomplex_t IT_0315 = 1.4142135623731*IT_0314;
    const ccomplex_t IT_0316 = 0.015625*IT_0315;
    const ccomplex_t IT_0317 = (-0.25)*IT_0084 + (-4)*IT_0086 + (-4)*IT_0093 +
       (-0.25)*IT_0095 + 0.25*IT_0097 + (-0.25)*IT_0100 + 0.25*IT_0104 + 0.25
      *IT_0109 + (-0.25)*IT_0111 + (-0.25)*IT_0136 + (-0.25)*IT_0149 + (-0.25)
      *IT_0159 + (-4)*IT_0161 + 0.25*IT_0168 + (-0.25)*IT_0170 + (-0.25)*IT_0172
       + (-0.25)*IT_0188 + 0.25*IT_0203 + 0.25*IT_0211 + (-0.25)*IT_0281;
    const ccomplex_t IT_0318 = IT_0080 + IT_0317;
    const ccomplex_t IT_0319 = IT_0045*IT_0318;
    const ccomplex_t IT_0320 = 1.4142135623731*IT_0319;
    const ccomplex_t IT_0321 = (-0.015625)*IT_0320;
    const ccomplex_t IT_0322 = IT_0032*IT_0040;
    const ccomplex_t IT_0323 = IT_0090 + IT_0102 + IT_0113 + IT_0130 + IT_0139
       + IT_0142 + IT_0158 + IT_0182 + IT_0198 + IT_0207 + IT_0229 + IT_0234 +
       IT_0236 + IT_0245 + IT_0248 + IT_0264 + IT_0267 + IT_0273 + IT_0276 +
       IT_0284 + IT_0287 + IT_0296 + IT_0302 + IT_0309 + IT_0311 + IT_0322;
    const ccomplex_t IT_0324 = IT_0032*IT_0038;
    const ccomplex_t IT_0325 = -IT_0107 + -IT_0121 + -IT_0124 + -IT_0133 + 16
      *IT_0146 + -IT_0152 + -IT_0156 + -IT_0157 + -IT_0166 + -IT_0192 + (-4)
      *IT_0219 + -IT_0223 + -IT_0226 + -IT_0232 + -IT_0242 + -IT_0250 + -IT_0254
       + -IT_0257 + -IT_0259 + -IT_0270 + -IT_0290 + -IT_0293 + -IT_0299 + 
      -IT_0305 + -IT_0307 + -IT_0324;
    const ccomplex_t IT_0326 = IT_0323 + IT_0325;
    const ccomplex_t IT_0327 = IT_0045*IT_0326;
    const ccomplex_t IT_0328 = 1.4142135623731*IT_0327;
    const ccomplex_t IT_0329 = 0.00390625*IT_0328;
    const ccomplex_t IT_0333 = (-0.25)*IT_0084 + (-4)*IT_0086 + (-4)*IT_0093 +
       (-0.25)*IT_0095 + (-0.25)*IT_0097 + 0.25*IT_0100 + (-0.25)*IT_0104 + (
      -0.25)*IT_0107 + 0.25*IT_0111 + (-0.25)*IT_0121 + 0.25*IT_0124 + (-0.25)
      *IT_0136 + (-0.25)*IT_0142 + (-0.25)*IT_0149 + (-0.25)*IT_0152 + (-0.25)
      *IT_0157 + (-0.25)*IT_0159 + (-4)*IT_0161 + (-0.25)*IT_0166 + 0.25*IT_0168
       + (-0.25)*IT_0170 + 0.25*IT_0172 + (-0.25)*IT_0182 + 0.25*IT_0192 + (
      -0.25)*IT_0207 + 0.25*IT_0211 + (-0.25)*IT_0264 + (-0.25)*IT_0281;
    const ccomplex_t IT_0334 = IT_0080 + IT_0333;
    const ccomplex_t IT_0335 = IT_0045*IT_0334;
    const ccomplex_t IT_0336 = 1.4142135623731*IT_0335;
    const ccomplex_t IT_0337 = (-0.015625)*IT_0336;
    const ccomplex_t IT_0338 = IT_0090 + IT_0109 + IT_0188 + IT_0198 + IT_0226
       + IT_0229 + IT_0234 + IT_0245 + IT_0259 + IT_0270 + IT_0276 + IT_0287 +
       IT_0290 + IT_0299 + IT_0302 + IT_0305 + IT_0307 + IT_0309 + IT_0322 +
       IT_0324;
    const ccomplex_t IT_0339 = -IT_0102 + -IT_0113 + -IT_0130 + -IT_0133 + 
      -IT_0139 + 16*IT_0146 + -IT_0156 + -IT_0158 + -IT_0203 + -IT_0208 + (-4)
      *IT_0219 + -IT_0223 + -IT_0232 + -IT_0236 + -IT_0242 + -IT_0248 + -IT_0250
       + -IT_0254 + -IT_0257 + -IT_0267 + -IT_0273 + -IT_0284 + -IT_0293 + 
      -IT_0296 + -IT_0311;
    const ccomplex_t IT_0340 = IT_0338 + IT_0339;
    const ccomplex_t IT_0341 = IT_0045*IT_0340;
    const ccomplex_t IT_0342 = 1.4142135623731*IT_0341;
    const ccomplex_t IT_0343 = 0.00390625*IT_0342;
    const ccomplex_t IT_0344 = IT_0071 + IT_0079;
    const ccomplex_t IT_0345 = (-0.25)*IT_0090 + (-0.25)*IT_0095 + 0.25
      *IT_0097 + 0.25*IT_0109 + 0.25*IT_0113 + 0.25*IT_0121 + 0.25*IT_0130 + (
      -0.25)*IT_0133 + (-0.25)*IT_0136 + 0.25*IT_0157 + (-0.25)*IT_0170 + (-0.25
      )*IT_0188 + (-0.25)*IT_0198 + (-0.25)*IT_0207 + (-0.25)*IT_0211 + (-0.25)
      *IT_0223 + (-0.25)*IT_0229 + 0.25*IT_0232 + (-0.25)*IT_0236 + (-0.25)
      *IT_0245 + (-0.25)*IT_0248 + (-0.25)*IT_0250 + 0.25*IT_0273 + (-0.25)
      *IT_0281 + 0.25*IT_0287 + (-0.25)*IT_0296 + (-0.25)*IT_0302 + 0.25*IT_0305
       + (-0.25)*IT_0307 + (-0.25)*IT_0324;
    const ccomplex_t IT_0346 = IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = IT_0045*IT_0346;
    const ccomplex_t IT_0348 = 1.4142135623731*IT_0347;
    const ccomplex_t IT_0349 = 0.015625*IT_0348;
    const ccomplex_t IT_0350 = IT_0084 + IT_0124 + IT_0142 + IT_0149 + IT_0156
       + IT_0159 + IT_0168 + IT_0182 + IT_0208 + IT_0226 + IT_0242 + IT_0257 +
       IT_0259 + IT_0264 + IT_0276 + IT_0284 + IT_0290 + IT_0293 + IT_0299 +
       IT_0309 + IT_0311 + IT_0322;
    const ccomplex_t IT_0351 = (-4)*IT_0064 + 16*IT_0086 + 16*IT_0093 + 
      -IT_0100 + -IT_0102 + -IT_0104 + -IT_0107 + -IT_0111 + -IT_0139 + 16
      *IT_0146 + -IT_0152 + -IT_0158 + 16*IT_0161 + -IT_0166 + -IT_0172 + 
      -IT_0192 + -IT_0203 + (-4)*IT_0219 + -IT_0234 + -IT_0254 + -IT_0267 + 
      -IT_0270;
    const ccomplex_t IT_0352 = IT_0350 + IT_0351;
    const ccomplex_t IT_0353 = IT_0045*IT_0352;
    const ccomplex_t IT_0354 = 1.4142135623731*IT_0353;
    const ccomplex_t IT_0355 = (-0.00390625)*IT_0354;
    const ccomplex_t IT_0356 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0004, IT_0007,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0357 = m_c*IT_0356;
    const ccomplex_t IT_0358 = (0 + _Complex_I*2)*IT_0004*IT_0009*IT_0010
      *IT_0017;
    const ccomplex_t IT_0359 = IT_0012*IT_0116*IT_0357*IT_0358;
    const ccomplex_t IT_0360 = 0.101321183642338*IT_0359;
    const ccomplex_t IT_0361 = mtylt_C0iC(9, 0, 0, 0, IT_0004, IT_0005,
       IT_0007, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0362 = 2*M_W*IT_0009*IT_0010;
    const ccomplex_t IT_0363 = -IT_0362;
    const ccomplex_t IT_0364 = IT_0361*IT_0363;
    const ccomplex_t IT_0365 = IT_0029*IT_0034*IT_0364;
    const ccomplex_t IT_0366 = 0.101321183642338*IT_0365;
    const ccomplex_t IT_0367 = mtylt_C0iC(9, 0, 0, 0, IT_0004, IT_0005,
       IT_0026, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0368 = IT_0363*IT_0367;
    const ccomplex_t IT_0369 = IT_0029*IT_0036*IT_0368;
    const ccomplex_t IT_0370 = 0.101321183642338*IT_0369;
    const ccomplex_t IT_0371 = 0.101321183642338*m_c*m_s;
    const ccomplex_t IT_0372 = mtylt_C0iC(0, 0, 0, 0, IT_0005, IT_0007,
       IT_0026, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0373 = IT_0020*IT_0021*IT_0060*IT_0372;
    const ccomplex_t IT_0374 = IT_0371*IT_0373;
    const ccomplex_t IT_0375 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0004, IT_0007,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0376 = -IT_0074;
    const ccomplex_t IT_0377 = IT_0375*IT_0376;
    const ccomplex_t IT_0378 = IT_0074*IT_0375;
    const ccomplex_t IT_0379 = (-8)*IT_0375;
    const ccomplex_t IT_0380 = Finite + IT_0379;
    const ccomplex_t IT_0381 = 2*IT_0074;
    const ccomplex_t IT_0382 = IT_0380*IT_0381;
    const ccomplex_t IT_0383 = (-2)*IT_0074;
    const ccomplex_t IT_0384 = IT_0380*IT_0383;
    const ccomplex_t IT_0385 = -IT_0378 + (-0.25)*IT_0382 + 0.25*IT_0384;
    const ccomplex_t IT_0386 = IT_0377 + IT_0385;
    const ccomplex_t IT_0387 = IT_0060*IT_0116*IT_0386;
    const ccomplex_t IT_0388 = 0.101321183642338*IT_0387;
    const ccomplex_t IT_0389 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0004, IT_0026,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0390 = IT_0376*IT_0389;
    const ccomplex_t IT_0391 = IT_0074*IT_0389;
    const ccomplex_t IT_0392 = (-8)*IT_0389;
    const ccomplex_t IT_0393 = Finite + IT_0392;
    const ccomplex_t IT_0394 = IT_0383*IT_0393;
    const ccomplex_t IT_0395 = IT_0381*IT_0393;
    const ccomplex_t IT_0396 = -IT_0391 + 0.25*IT_0394 + (-0.25)*IT_0395;
    const ccomplex_t IT_0397 = IT_0390 + IT_0396;
    const ccomplex_t IT_0398 = IT_0060*IT_0075*IT_0397;
    const ccomplex_t IT_0399 = 0.101321183642338*IT_0398;
    const ccomplex_t IT_0400 = mtylt_C0iC(0, 0, 0, 0, IT_0004, IT_0007,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0401 = m_c*IT_0400;
    const ccomplex_t IT_0402 = (0 + _Complex_I*2)*IT_0004*IT_0009*IT_0010;
    const ccomplex_t IT_0403 = IT_0060*IT_0177*IT_0401*IT_0402;
    const ccomplex_t IT_0404 = 0.101321183642338*IT_0403;
    const ccomplex_t IT_0405 = mtylt_C0iC(0, 0, 0, 0, IT_0004, IT_0173,
       IT_0026, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0406 = m_s*IT_0405;
    const ccomplex_t IT_0407 = IT_0060*IT_0239*IT_0402*IT_0406;
    const ccomplex_t IT_0408 = 0.101321183642338*IT_0407;
    const ccomplex_t IT_0409 = mtylt_C0iC(9, 0, 0, 0, IT_0004, IT_0007,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0410 = (0 + _Complex_I*2)*M_W*IT_0009*IT_0010;
    const ccomplex_t IT_0411 = -IT_0410;
    const ccomplex_t IT_0412 = IT_0409*IT_0411;
    const ccomplex_t IT_0413 = IT_0012*IT_0177*IT_0412;
    const ccomplex_t IT_0414 = 0.101321183642338*IT_0413;
    const ccomplex_t IT_0415 = mtylt_C0iC(9, 0, 0, 0, IT_0004, IT_0173,
       IT_0026, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0416 = IT_0411*IT_0415;
    const ccomplex_t IT_0417 = IT_0012*IT_0239*IT_0416;
    const ccomplex_t IT_0418 = 0.101321183642338*IT_0417;
    const ccomplex_t IT_0419 = 0.101321183642338*m_s;
    const ccomplex_t IT_0420 = mtylt_C0iC(0, 0, 0, 0, IT_0004, IT_0005,
       IT_0026, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0421 = (0 + _Complex_I*2)*IT_0004*IT_0009*IT_0010
      *IT_0016*IT_0018;
    const ccomplex_t IT_0422 = -IT_0421;
    const ccomplex_t IT_0423 = IT_0012*IT_0055*IT_0420*IT_0422;
    const ccomplex_t IT_0424 = IT_0419*IT_0423;
    const ccomplex_t IT_0425 = 0.101321183642338*m_c;
    const ccomplex_t IT_0426 = mtylt_C0iC(0, 0, 0, 0, IT_0004, IT_0005,
       IT_0007, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0427 = IT_0029*IT_0066*IT_0422*IT_0426;
    const ccomplex_t IT_0428 = IT_0425*IT_0427;
    const ccomplex_t IT_0429 = mtylt_C0iC(9, 0, 0, 0, IT_0005, IT_0007,
       IT_0026, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0430 = (-2)*IT_0429;
    const ccomplex_t IT_0431 = Finite + IT_0430;
    const ccomplex_t IT_0432 = IT_0055*IT_0060*IT_0066*IT_0431;
    const ccomplex_t IT_0433 = 0.101321183642338*IT_0432;
    const ccomplex_t IT_0434 = IT_0029*IT_0116*IT_0357*IT_0358;
    const ccomplex_t IT_0435 = 0.101321183642338*IT_0434;
    const ccomplex_t IT_0436 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0004, IT_0026,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0437 = m_s*IT_0436;
    const ccomplex_t IT_0438 = IT_0012*IT_0075*IT_0358*IT_0437;
    const ccomplex_t IT_0439 = 0.101321183642338*IT_0438;
    const ccomplex_t IT_0440 = IT_0029*IT_0075*IT_0358*IT_0437;
    const ccomplex_t IT_0441 = 0.101321183642338*IT_0440;
    const ccomplex_t IT_0442 = IT_0029*IT_0177*IT_0412;
    const ccomplex_t IT_0443 = 0.101321183642338*IT_0442;
    const ccomplex_t IT_0444 = IT_0029*IT_0239*IT_0416;
    const ccomplex_t IT_0445 = 0.101321183642338*IT_0444;
    const ccomplex_t IT_0446 = IT_0012*IT_0020*IT_0422*IT_0426;
    const ccomplex_t IT_0447 = IT_0425*IT_0446;
    const ccomplex_t IT_0448 = IT_0021*IT_0029*IT_0420*IT_0422;
    const ccomplex_t IT_0449 = IT_0419*IT_0448;
    const ccomplex_t IT_0450 = IT_0012*IT_0034*IT_0364;
    const ccomplex_t IT_0451 = 0.101321183642338*IT_0450;
    const ccomplex_t IT_0452 = IT_0012*IT_0036*IT_0368;
    const ccomplex_t IT_0453 = 0.101321183642338*IT_0452;
    const ccomplex_t IT_0454 = (0 + _Complex_I*2)*M_W*IT_0009*IT_0010*IT_0015;
    const ccomplex_t IT_0455 = -IT_0454;
    const ccomplex_t IT_0456 = IT_0361*IT_0455;
    const ccomplex_t IT_0457 = IT_0361*IT_0454;
    const ccomplex_t IT_0458 = (-8)*IT_0361;
    const ccomplex_t IT_0459 = Finite + IT_0458;
    const ccomplex_t IT_0460 = IT_0454*IT_0459;
    const ccomplex_t IT_0461 = IT_0455*IT_0459;
    const ccomplex_t IT_0462 = -IT_0457 + (-0.5)*IT_0460 + 0.5*IT_0461;
    const ccomplex_t IT_0463 = IT_0456 + IT_0462;
    const ccomplex_t IT_0464 = IT_0060*IT_0066*IT_0463;
    const ccomplex_t IT_0465 = 0.101321183642338*IT_0464;
    const ccomplex_t IT_0466 = IT_0367*IT_0455;
    const ccomplex_t IT_0467 = (-8)*IT_0367;
    const ccomplex_t IT_0468 = Finite + IT_0467;
    const ccomplex_t IT_0469 = IT_0454*IT_0468;
    const ccomplex_t IT_0470 = IT_0455*IT_0468;
    const ccomplex_t IT_0471 = IT_0367*IT_0454;
    const ccomplex_t IT_0472 = (-0.5)*IT_0469 + 0.5*IT_0470 + -IT_0471;
    const ccomplex_t IT_0473 = IT_0466 + IT_0472;
    const ccomplex_t IT_0474 = IT_0055*IT_0060*IT_0473;
    const ccomplex_t IT_0475 = 0.101321183642338*IT_0474;
    const ccomplex_t IT_0476 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0007, IT_0026,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0477 = (-8)*IT_0476;
    const ccomplex_t IT_0478 = Finite + IT_0477;
    const ccomplex_t IT_0479 = (-4)*IT_0476;
    const ccomplex_t IT_0480 = Finite + IT_0479;
    const ccomplex_t IT_0481 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0007, IT_0026,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0482 = IT_0126*IT_0481;
    const ccomplex_t IT_0483 = 0.5*IT_0478 + -IT_0480 + (-0.5)*IT_0482;
    const ccomplex_t IT_0484 = IT_0476 + IT_0483;
    const ccomplex_t IT_0485 = IT_0060*IT_0075*IT_0116;
    const ccomplex_t IT_0486 = 0.101321183642338*IT_0485;
    const ccomplex_t IT_0487 = IT_0484*IT_0486;
    const ccomplex_t IT_0488 = mtylt_C0iC(0, 0, 0, 0, IT_0007, IT_0173,
       IT_0026, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0489 = IT_0126*IT_0488;
    const ccomplex_t IT_0490 = mtylt_C0iC(9, 0, 0, 0, IT_0007, IT_0173,
       IT_0026, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0491 = (-8)*IT_0490;
    const ccomplex_t IT_0492 = Finite + IT_0491;
    const ccomplex_t IT_0493 = 2*IT_0490 + 0.5*IT_0492;
    const ccomplex_t IT_0494 = IT_0489 + IT_0493;
    const ccomplex_t IT_0495 = IT_0060*IT_0177*IT_0239;
    const ccomplex_t IT_0496 = 0.101321183642338*IT_0495;
    const ccomplex_t IT_0497 = IT_0494*IT_0496;
    const ccomplex_t IT_0498 = (-8)*IT_0429;
    const ccomplex_t IT_0499 = Finite + IT_0498;
    const ccomplex_t IT_0500 = IT_0126*IT_0372;
    const ccomplex_t IT_0501 = 0.25*IT_0499 + 0.5*IT_0500;
    const ccomplex_t IT_0502 = IT_0429 + IT_0501;
    const ccomplex_t IT_0503 = IT_0034*IT_0036*IT_0060;
    const ccomplex_t IT_0504 = 0.101321183642338*IT_0503;
    const ccomplex_t IT_0505 = IT_0502*IT_0504;
    const ccomplex_t IT_0506 = IT_0360 + -IT_0366 + IT_0370 + (-2)*IT_0374 + 2
      *IT_0388 + (-2)*IT_0399 + IT_0404 + IT_0408 + IT_0414 + -IT_0418 + IT_0424
       + IT_0428 + (-2)*IT_0433 + IT_0435 + IT_0439 + IT_0441 + IT_0443 + 
      -IT_0445 + IT_0447 + IT_0449 + IT_0451 + -IT_0453 + 2*IT_0465 + (-2)
      *IT_0475 + 4*IT_0487 + -IT_0497 + 2*IT_0505;
    const ccomplex_t IT_0507 = cpow(2*s_13 + IT_0004 + -IT_0006 + -IT_0007 + 
      -reg_prop, -1);
    const ccomplex_t IT_0508 = IT_0506*IT_0507;
    const ccomplex_t IT_0509 = IT_0045*IT_0062*IT_0508;
    const ccomplex_t IT_0510 = 1.4142135623731*IT_0509;
    const ccomplex_t IT_0511 = (-0.00390625)*IT_0510;
    const ccomplex_t IT_0512 = -IT_0508;
    const ccomplex_t IT_0513 = -IT_0512;
    const ccomplex_t IT_0514 = IT_0045*IT_0062*IT_0513;
    const ccomplex_t IT_0515 = 1.4142135623731*IT_0514;
    const ccomplex_t IT_0516 = 0.00390625*IT_0515;
    const ccomplex_t IT_0517 = 0.5*IT_0478 + -IT_0480 + 0.5*IT_0482;
    const ccomplex_t IT_0518 = IT_0476 + IT_0517;
    const ccomplex_t IT_0519 = IT_0486*IT_0518;
    const ccomplex_t IT_0520 = 0.25*IT_0499 + (-0.5)*IT_0500;
    const ccomplex_t IT_0521 = IT_0429 + IT_0520;
    const ccomplex_t IT_0522 = IT_0504*IT_0521;
    const ccomplex_t IT_0523 = (-2)*IT_0490 + (-0.5)*IT_0492;
    const ccomplex_t IT_0524 = IT_0489 + IT_0523;
    const ccomplex_t IT_0525 = IT_0496*IT_0524;
    const ccomplex_t IT_0526 = IT_0360 + IT_0366 + IT_0370 + (-2)*IT_0374 + (
      -2)*IT_0388 + 2*IT_0399 + -IT_0404 + -IT_0408 + IT_0414 + IT_0418 + 
      -IT_0424 + -IT_0428 + 2*IT_0433 + -IT_0435 + -IT_0439 + IT_0441 + -IT_0443
       + -IT_0445 + IT_0447 + IT_0449 + IT_0451 + IT_0453 + (-2)*IT_0465 + 2
      *IT_0475 + (-4)*IT_0519 + 2*IT_0522 + IT_0525;
    const ccomplex_t IT_0527 = IT_0507*IT_0526;
    const ccomplex_t IT_0528 = -IT_0527;
    const ccomplex_t IT_0529 = -IT_0528;
    const ccomplex_t IT_0530 = IT_0045*IT_0062*IT_0529;
    const ccomplex_t IT_0531 = 1.4142135623731*IT_0530;
    const ccomplex_t IT_0532 = (-0.00390625)*IT_0531;
    const ccomplex_t IT_0533 = 0.00390625*IT_0531;
    const ccomplex_t IT_0534 = IT_0014*IT_0116*IT_0357*IT_0358;
    const ccomplex_t IT_0535 = 0.101321183642338*IT_0534;
    const ccomplex_t IT_0536 = IT_0042*IT_0116*IT_0357*IT_0358;
    const ccomplex_t IT_0537 = 0.101321183642338*IT_0536;
    const ccomplex_t IT_0538 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0004, IT_0006,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0539 = m_b*IT_0538;
    const ccomplex_t IT_0540 = IT_0042*IT_0075*IT_0358*IT_0539;
    const ccomplex_t IT_0541 = 0.101321183642338*IT_0540;
    const ccomplex_t IT_0542 = mtylt_C0iC(0, 0, 0, 0, IT_0004, IT_0006,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0543 = m_b*IT_0542;
    const ccomplex_t IT_0544 = IT_0062*IT_0179*IT_0402*IT_0543;
    const ccomplex_t IT_0545 = 0.101321183642338*IT_0544;
    const ccomplex_t IT_0546 = IT_0062*IT_0177*IT_0401*IT_0402;
    const ccomplex_t IT_0547 = 0.101321183642338*IT_0546;
    const ccomplex_t IT_0548 = mtylt_C0iC(9, 0, 0, 0, IT_0004, IT_0006,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0549 = IT_0411*IT_0548;
    const ccomplex_t IT_0550 = IT_0014*IT_0179*IT_0549;
    const ccomplex_t IT_0551 = 0.101321183642338*IT_0550;
    const ccomplex_t IT_0552 = IT_0014*IT_0177*IT_0412;
    const ccomplex_t IT_0553 = 0.101321183642338*IT_0552;
    const ccomplex_t IT_0554 = IT_0042*IT_0179*IT_0549;
    const ccomplex_t IT_0555 = 0.101321183642338*IT_0554;
    const ccomplex_t IT_0556 = IT_0042*IT_0177*IT_0412;
    const ccomplex_t IT_0557 = 0.101321183642338*IT_0556;
    const ccomplex_t IT_0558 = IT_0014*IT_0066*IT_0422*IT_0426;
    const ccomplex_t IT_0559 = IT_0425*IT_0558;
    const ccomplex_t IT_0560 = 0.101321183642338*m_b;
    const ccomplex_t IT_0561 = mtylt_C0iC(0, 0, 0, 0, IT_0004, IT_0005,
       IT_0006, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0562 = IT_0014*IT_0021*IT_0422*IT_0561;
    const ccomplex_t IT_0563 = IT_0560*IT_0562;
    const ccomplex_t IT_0564 = IT_0042*IT_0055*IT_0422*IT_0561;
    const ccomplex_t IT_0565 = IT_0560*IT_0564;
    const ccomplex_t IT_0566 = IT_0020*IT_0042*IT_0422*IT_0426;
    const ccomplex_t IT_0567 = IT_0425*IT_0566;
    const ccomplex_t IT_0568 = mtylt_C0iC(9, 0, 0, 0, IT_0004, IT_0005,
       IT_0006, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0569 = IT_0362*IT_0568;
    const ccomplex_t IT_0570 = IT_0014*IT_0185*IT_0569;
    const ccomplex_t IT_0571 = 0.101321183642338*IT_0570;
    const ccomplex_t IT_0572 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0006, IT_0007,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0573 = (-8)*IT_0572;
    const ccomplex_t IT_0574 = Finite + IT_0573;
    const ccomplex_t IT_0575 = (-4)*IT_0572;
    const ccomplex_t IT_0576 = Finite + IT_0575;
    const ccomplex_t IT_0577 = mtylt_C0iC(0, 0, 0, 0, 0, IT_0006, IT_0007,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0578 = IT_0154*IT_0577;
    const ccomplex_t IT_0579 = 0.5*IT_0574 + -IT_0576 + (-0.5)*IT_0578;
    const ccomplex_t IT_0580 = IT_0572 + IT_0579;
    const ccomplex_t IT_0581 = IT_0062*IT_0075*IT_0116;
    const ccomplex_t IT_0582 = 0.101321183642338*IT_0581;
    const ccomplex_t IT_0583 = IT_0580*IT_0582;
    const ccomplex_t IT_0584 = mtylt_C0iC(9, 0, 0, 0, IT_0005, IT_0006,
       IT_0007, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0585 = (-8)*IT_0584;
    const ccomplex_t IT_0586 = Finite + IT_0585;
    const ccomplex_t IT_0587 = mtylt_C0iC(0, 0, 0, 0, IT_0005, IT_0006,
       IT_0007, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0588 = IT_0154*IT_0587;
    const ccomplex_t IT_0589 = 0.25*IT_0586 + 0.5*IT_0588;
    const ccomplex_t IT_0590 = IT_0584 + IT_0589;
    const ccomplex_t IT_0591 = IT_0034*IT_0062*IT_0185;
    const ccomplex_t IT_0592 = 0.101321183642338*IT_0591;
    const ccomplex_t IT_0593 = IT_0590*IT_0592;
    const ccomplex_t IT_0594 = IT_0042*IT_0185*IT_0569;
    const ccomplex_t IT_0595 = 0.101321183642338*IT_0594;
    const ccomplex_t IT_0596 = IT_0014*IT_0075*IT_0358*IT_0539;
    const ccomplex_t IT_0597 = 0.101321183642338*IT_0596;
    const ccomplex_t IT_0598 = mtylt_C0iC(0, 0, 0, 0, IT_0006, IT_0007,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0599 = IT_0154*IT_0598;
    const ccomplex_t IT_0600 = mtylt_C0iC(9, 0, 0, 0, IT_0006, IT_0007,
       IT_0173, getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0601 = (-8)*IT_0600;
    const ccomplex_t IT_0602 = Finite + IT_0601;
    const ccomplex_t IT_0603 = 2*IT_0600 + 0.5*IT_0602;
    const ccomplex_t IT_0604 = IT_0599 + IT_0603;
    const ccomplex_t IT_0605 = IT_0062*IT_0177*IT_0179;
    const ccomplex_t IT_0606 = 0.101321183642338*IT_0605;
    const ccomplex_t IT_0607 = IT_0604*IT_0606;
    const ccomplex_t IT_0608 = IT_0062*IT_0066*IT_0463;
    const ccomplex_t IT_0609 = 0.101321183642338*IT_0608;
    const ccomplex_t IT_0610 = IT_0455*IT_0568;
    const ccomplex_t IT_0611 = IT_0454*IT_0568;
    const ccomplex_t IT_0612 = (-8)*IT_0568;
    const ccomplex_t IT_0613 = Finite + IT_0612;
    const ccomplex_t IT_0614 = IT_0454*IT_0613;
    const ccomplex_t IT_0615 = IT_0455*IT_0613;
    const ccomplex_t IT_0616 = -IT_0611 + (-0.5)*IT_0614 + 0.5*IT_0615;
    const ccomplex_t IT_0617 = IT_0610 + IT_0616;
    const ccomplex_t IT_0618 = IT_0055*IT_0062*IT_0617;
    const ccomplex_t IT_0619 = 0.101321183642338*IT_0618;
    const ccomplex_t IT_0620 = IT_0375*IT_0381;
    const ccomplex_t IT_0621 = IT_0376*IT_0380;
    const ccomplex_t IT_0622 = 0.25*IT_0382 + (-0.5)*IT_0621;
    const ccomplex_t IT_0623 = IT_0620 + IT_0622;
    const ccomplex_t IT_0624 = IT_0062*IT_0116*IT_0623;
    const ccomplex_t IT_0625 = 0.101321183642338*IT_0624;
    const ccomplex_t IT_0626 = IT_0361*IT_0362;
    const ccomplex_t IT_0627 = IT_0014*IT_0034*IT_0626;
    const ccomplex_t IT_0628 = 0.101321183642338*IT_0627;
    const ccomplex_t IT_0629 = (-2)*IT_0584;
    const ccomplex_t IT_0630 = Finite + IT_0629;
    const ccomplex_t IT_0631 = IT_0055*IT_0062*IT_0066*IT_0630;
    const ccomplex_t IT_0632 = 0.101321183642338*IT_0631;
    const ccomplex_t IT_0633 = 0.101321183642338*m_b*m_c;
    const ccomplex_t IT_0634 = IT_0020*IT_0021*IT_0062*IT_0587;
    const ccomplex_t IT_0635 = IT_0633*IT_0634;
    const ccomplex_t IT_0636 = mtylt_C0iC(9, 0, 0, 0, 0, IT_0004, IT_0006,
       getIntegrationParameters()->reg_int);
    const ccomplex_t IT_0637 = (-8)*IT_0636;
    const ccomplex_t IT_0638 = Finite + IT_0637;
    const ccomplex_t IT_0639 = IT_0376*IT_0638;
    const ccomplex_t IT_0640 = IT_0381*IT_0638;
    const ccomplex_t IT_0641 = (-0.5)*IT_0639 + 0.25*IT_0640;
    const ccomplex_t IT_0642 = IT_0381*IT_0636;
    const ccomplex_t IT_0643 = IT_0641 + IT_0642;
    const ccomplex_t IT_0644 = IT_0062*IT_0075*IT_0643;
    const ccomplex_t IT_0645 = 0.101321183642338*IT_0644;
    const ccomplex_t IT_0646 = IT_0034*IT_0042*IT_0626;
    const ccomplex_t IT_0647 = 0.101321183642338*IT_0646;
    const ccomplex_t IT_0648 = IT_0535 + IT_0537 + IT_0541 + IT_0545 + IT_0547
       + -IT_0551 + IT_0553 + -IT_0555 + IT_0557 + IT_0559 + IT_0563 + IT_0565 +
       IT_0567 + -IT_0571 + 4*IT_0583 + 2*IT_0593 + IT_0595 + IT_0597 + -IT_0607
       + 2*IT_0609 + (-2)*IT_0619 + (-2)*IT_0625 + IT_0628 + (-2)*IT_0632 + (-2)
      *IT_0635 + 2*IT_0645 + -IT_0647;
    const ccomplex_t IT_0649 = cpow((-2)*s_24 + IT_0004 + -IT_0007 + -IT_0026 
      + -reg_prop, -1);
    const ccomplex_t IT_0650 = IT_0648*IT_0649;
    const ccomplex_t IT_0651 = IT_0045*IT_0060*IT_0650;
    const ccomplex_t IT_0652 = 1.4142135623731*IT_0651;
    const ccomplex_t IT_0653 = (-0.00390625)*IT_0652;
    const ccomplex_t IT_0654 = 0.5*IT_0574 + -IT_0576 + 0.5*IT_0578;
    const ccomplex_t IT_0655 = IT_0572 + IT_0654;
    const ccomplex_t IT_0656 = IT_0582*IT_0655;
    const ccomplex_t IT_0657 = (-2)*IT_0600 + (-0.5)*IT_0602;
    const ccomplex_t IT_0658 = IT_0599 + IT_0657;
    const ccomplex_t IT_0659 = IT_0606*IT_0658;
    const ccomplex_t IT_0660 = 0.25*IT_0586 + (-0.5)*IT_0588;
    const ccomplex_t IT_0661 = IT_0584 + IT_0660;
    const ccomplex_t IT_0662 = IT_0592*IT_0661;
    const ccomplex_t IT_0663 = IT_0535 + -IT_0537 + IT_0541 + IT_0545 +
       IT_0547 + IT_0551 + IT_0553 + -IT_0555 + -IT_0557 + IT_0559 + -IT_0563 +
       IT_0565 + -IT_0567 + IT_0571 + IT_0595 + -IT_0597 + 2*IT_0609 + (-2)
      *IT_0619 + (-2)*IT_0625 + IT_0628 + (-2)*IT_0632 + 2*IT_0635 + 2*IT_0645 +
       IT_0647 + 4*IT_0656 + -IT_0659 + (-2)*IT_0662;
    const ccomplex_t IT_0664 = IT_0649*IT_0663;
    const ccomplex_t IT_0665 = (-0.5)*IT_0664;
    const ccomplex_t IT_0666 = (-2)*IT_0665;
    const ccomplex_t IT_0667 = IT_0045*IT_0060*IT_0666;
    const ccomplex_t IT_0668 = 1.4142135623731*IT_0667;
    const ccomplex_t IT_0669 = 0.00390625*IT_0668;
    const ccomplex_t IT_0670 = 0.00390625*IT_0652;
    const ccomplex_t IT_0671 = -IT_0664;
    const ccomplex_t IT_0672 = -IT_0671;
    const ccomplex_t IT_0673 = IT_0045*IT_0060*IT_0672;
    const ccomplex_t IT_0674 = 1.4142135623731*IT_0673;
    const ccomplex_t IT_0675 = (-0.00390625)*IT_0674;
    const ccomplex_t IT_0676 = (0 + _Complex_I*0.0055242717280199)*IT_0000
      *IT_0001*IT_0002*IT_0032*IT_0038 + (0 + _Complex_I*(-0.0055242717280199))
      *IT_0000*IT_0001*IT_0002*IT_0032*IT_0040 + (0 + _Complex_I*(
      -0.0055242717280199))*IT_0000*IT_0001*IT_0002*IT_0032*IT_0044 + (0 +
       _Complex_I*1)*IT_0216 + (0 + _Complex_I*1)*IT_0316 + (0 + _Complex_I*-1)
      *IT_0321 + (0 + _Complex_I*-1)*IT_0329 + (0 + _Complex_I*-1)*IT_0337 + (0 
      + _Complex_I*-1)*IT_0343 + (0 + _Complex_I*1)*IT_0349 + (0 + _Complex_I*1)
      *IT_0355 + IT_0511 + -IT_0516 + -IT_0532 + IT_0533 + IT_0653 + -IT_0669 + 
      -IT_0670 + IT_0675;
    return create_ccomplex_return(IT_0676);
}

