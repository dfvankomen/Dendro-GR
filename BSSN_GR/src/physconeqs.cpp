// Dendro: {{{
// Dendro: original ops: 1298272
// Dendro: printing temp variables
const double DENDRO_0  = grad_0_Gt0[pp];
const double DENDRO_1  = grad_0_Gt1[pp];
const double DENDRO_2  = 4 * gt1[pp];
const double DENDRO_3  = grad_0_Gt2[pp];
const double DENDRO_4  = 4 * gt2[pp];
const double DENDRO_5  = pow(chi[pp], -2);
const double DENDRO_6  = grad_0_chi[pp];
const double DENDRO_7  = (DENDRO_6 * DENDRO_6);
const double DENDRO_8  = DENDRO_5 * DENDRO_7;
const double DENDRO_9  = gt1[pp] * gt2[pp];
const double DENDRO_10 = -DENDRO_9 + gt0[pp] * gt4[pp];
const double DENDRO_11 = pow(gt4[pp], 2);
const double DENDRO_12 = pow(gt1[pp], 2);
const double DENDRO_13 = pow(gt2[pp], 2);
const double DENDRO_14 = gt0[pp] * gt3[pp];
const double DENDRO_15 = DENDRO_11 * gt0[pp] + DENDRO_12 * gt5[pp] +
                         DENDRO_13 * gt3[pp] - DENDRO_14 * gt5[pp] -
                         2 * DENDRO_9 * gt4[pp];
const double DENDRO_16 = 1.0 / DENDRO_15;
const double DENDRO_17 = DENDRO_10 * DENDRO_16;
const double DENDRO_18 = 4.0 * DENDRO_17;
const double DENDRO_19 = DENDRO_18 * grad2_1_2_gt0[pp];
const double DENDRO_20 = gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp];
const double DENDRO_21 = DENDRO_16 * DENDRO_20;
const double DENDRO_22 = 4.0 * DENDRO_21;
const double DENDRO_23 = DENDRO_22 * grad2_0_1_gt0[pp];
const double DENDRO_24 = grad2_0_2_gt0[pp];
const double DENDRO_25 = gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp];
const double DENDRO_26 = DENDRO_16 * DENDRO_25;
const double DENDRO_27 = 4 * DENDRO_26;
const double DENDRO_28 = grad2_2_2_gt0[pp];
const double DENDRO_29 = -DENDRO_12 + DENDRO_14;
const double DENDRO_30 = DENDRO_16 * DENDRO_29;
const double DENDRO_31 = 2.0 * DENDRO_30;
const double DENDRO_32 = grad2_1_1_gt0[pp];
const double DENDRO_33 = -DENDRO_13 + gt0[pp] * gt5[pp];
const double DENDRO_34 = DENDRO_16 * DENDRO_33;
const double DENDRO_35 = 2.0 * DENDRO_34;
const double DENDRO_36 = grad2_0_0_gt0[pp];
const double DENDRO_37 = -DENDRO_11 + gt3[pp] * gt5[pp];
const double DENDRO_38 = DENDRO_16 * DENDRO_37;
const double DENDRO_39 = 2.0 * DENDRO_38;
const double DENDRO_40 = grad_2_gt0[pp];
const double DENDRO_41 = grad_0_gt5[pp];
const double DENDRO_42 = grad_0_gt4[pp];
const double DENDRO_43 = grad_2_gt1[pp];
const double DENDRO_44 = grad_1_gt2[pp];
const double DENDRO_45 = DENDRO_42 + DENDRO_43 - DENDRO_44;
const double DENDRO_46 =
    -DENDRO_20 * DENDRO_45 + DENDRO_25 * DENDRO_41 + DENDRO_37 * DENDRO_40;
const double DENDRO_47 = -DENDRO_46;
const double DENDRO_48 = DENDRO_40 * DENDRO_47;
const double DENDRO_49 = pow(DENDRO_15, -2);
const double DENDRO_50 = DENDRO_29 * DENDRO_49;
const double DENDRO_51 = 3.0 * DENDRO_50;
const double DENDRO_52 = grad_1_gt0[pp];
const double DENDRO_53 = grad_0_gt3[pp];
const double DENDRO_54 = DENDRO_42 - DENDRO_43 + DENDRO_44;
const double DENDRO_55 =
    -DENDRO_20 * DENDRO_53 + DENDRO_25 * DENDRO_54 + DENDRO_37 * DENDRO_52;
const double DENDRO_56 = -DENDRO_55;
const double DENDRO_57 = DENDRO_52 * DENDRO_56;
const double DENDRO_58 = DENDRO_33 * DENDRO_49;
const double DENDRO_59 = 3.0 * DENDRO_58;
const double DENDRO_60 = grad_0_gt0[pp];
const double DENDRO_61 = 0.5 * DENDRO_60;
const double DENDRO_62 = grad_0_gt2[pp];
const double DENDRO_63 = 1.0 * DENDRO_62;
const double DENDRO_64 = 0.5 * DENDRO_40;
const double DENDRO_65 = DENDRO_63 - DENDRO_64;
const double DENDRO_66 = grad_0_gt1[pp];
const double DENDRO_67 = 1.0 * DENDRO_66;
const double DENDRO_68 = 0.5 * DENDRO_52;
const double DENDRO_69 = DENDRO_67 - DENDRO_68;
const double DENDRO_70 =
    -DENDRO_20 * DENDRO_69 + DENDRO_25 * DENDRO_65 + DENDRO_37 * DENDRO_61;
const double DENDRO_71 = -DENDRO_70;
const double DENDRO_72 = DENDRO_37 * DENDRO_71;
const double DENDRO_73 = 6.0 * DENDRO_49;
const double DENDRO_74 = 0.25 * DENDRO_53;
const double DENDRO_75 = grad_1_gt1[pp];
const double DENDRO_76 = 1.0 * DENDRO_75;
const double DENDRO_77 = -DENDRO_76;
const double DENDRO_78 = DENDRO_74 + DENDRO_77;
const double DENDRO_79 = DENDRO_10 * DENDRO_54 + DENDRO_20 * DENDRO_52;
const double DENDRO_80 = -DENDRO_33 * DENDRO_53 + DENDRO_79;
const double DENDRO_81 = 4 * DENDRO_58;
const double DENDRO_82 = 0.25 * DENDRO_41;
const double DENDRO_83 = grad_2_gt2[pp];
const double DENDRO_84 = 1.0 * DENDRO_83;
const double DENDRO_85 = -DENDRO_84;
const double DENDRO_86 = DENDRO_82 + DENDRO_85;
const double DENDRO_87 =
    -DENDRO_10 * DENDRO_45 + DENDRO_25 * DENDRO_40 + DENDRO_29 * DENDRO_41;
const double DENDRO_88 = -DENDRO_87;
const double DENDRO_89 = 4 * DENDRO_50;
const double DENDRO_90 = DENDRO_10 * DENDRO_41 + DENDRO_20 * DENDRO_40;
const double DENDRO_91 = -DENDRO_33 * DENDRO_45 + DENDRO_90;
const double DENDRO_92 = 0.25 * DENDRO_42;
const double DENDRO_93 = -DENDRO_92;
const double DENDRO_94 = 0.25 * DENDRO_44;
const double DENDRO_95 = 0.75 * DENDRO_43;
const double DENDRO_96 =
    DENDRO_89 * DENDRO_91 * (DENDRO_93 + DENDRO_94 + DENDRO_95);
const double DENDRO_97 =
    -DENDRO_10 * DENDRO_53 + DENDRO_25 * DENDRO_52 + DENDRO_29 * DENDRO_54;
const double DENDRO_98  = -DENDRO_97;
const double DENDRO_99  = 0.75 * DENDRO_44;
const double DENDRO_100 = 0.25 * DENDRO_43;
const double DENDRO_101 = DENDRO_100 + DENDRO_93 + DENDRO_99;
const double DENDRO_102 = DENDRO_10 * DENDRO_65 + DENDRO_20 * DENDRO_61;
const double DENDRO_103 = DENDRO_102 - DENDRO_33 * DENDRO_69;
const double DENDRO_104 = DENDRO_103 * DENDRO_37;
const double DENDRO_105 = 4 * DENDRO_49;
const double DENDRO_106 = DENDRO_104 * DENDRO_105 * (DENDRO_67 + DENDRO_68);
const double DENDRO_107 = DENDRO_63 + DENDRO_64;
const double DENDRO_108 =
    -DENDRO_10 * DENDRO_69 + DENDRO_25 * DENDRO_61 + DENDRO_29 * DENDRO_65;
const double DENDRO_109 = -DENDRO_108;
const double DENDRO_110 = DENDRO_109 * DENDRO_37;
const double DENDRO_111 = 0.25 * DENDRO_40;
const double DENDRO_112 = DENDRO_111 * DENDRO_56;
const double DENDRO_113 = DENDRO_10 * DENDRO_49;
const double DENDRO_114 = 4 * DENDRO_113;
const double DENDRO_115 = 0.25 * DENDRO_52;
const double DENDRO_116 = DENDRO_115 * DENDRO_47;
const double DENDRO_117 = DENDRO_52 * DENDRO_80;
const double DENDRO_118 = DENDRO_103 * DENDRO_53;
const double DENDRO_119 = DENDRO_117 + DENDRO_118;
const double DENDRO_120 = DENDRO_20 * DENDRO_49;
const double DENDRO_121 = 2.0 * DENDRO_120;
const double DENDRO_122 = DENDRO_109 * DENDRO_41;
const double DENDRO_123 = DENDRO_25 * DENDRO_49;
const double DENDRO_124 = 2.0 * DENDRO_123;
const double DENDRO_125 = DENDRO_47 * DENDRO_60;
const double DENDRO_126 = DENDRO_40 * DENDRO_71;
const double DENDRO_127 = DENDRO_56 * DENDRO_60;
const double DENDRO_128 = DENDRO_52 * DENDRO_71;
const double DENDRO_129 = 0.25 * DENDRO_125;
const double DENDRO_130 = 4 * DENDRO_123;
const double DENDRO_131 = 0.25 * DENDRO_127;
const double DENDRO_132 = 4 * DENDRO_120;
const double DENDRO_133 = DENDRO_53 * DENDRO_91;
const double DENDRO_134 = 0.25 * DENDRO_133;
const double DENDRO_135 = -DENDRO_42 + DENDRO_43 + DENDRO_44;
const double DENDRO_136 = 0.5 * DENDRO_135;
const double DENDRO_137 = DENDRO_134 + DENDRO_136 * DENDRO_80;
const double DENDRO_138 = DENDRO_82 * DENDRO_98;
const double DENDRO_139 = DENDRO_52 * DENDRO_91;
const double DENDRO_140 = DENDRO_103 * DENDRO_45;
const double DENDRO_141 = DENDRO_124 * (DENDRO_139 + DENDRO_140);
const double DENDRO_142 = DENDRO_109 * DENDRO_54;
const double DENDRO_143 = 0.5 * DENDRO_53;
const double DENDRO_144 = DENDRO_143 + DENDRO_77;
const double DENDRO_145 = DENDRO_144 * DENDRO_91;
const double DENDRO_146 = DENDRO_45 * DENDRO_80;
const double DENDRO_147 = 0.25 * DENDRO_146;
const double DENDRO_148 = -DENDRO_147;
const double DENDRO_149 = 0.5 * DENDRO_41;
const double DENDRO_150 = DENDRO_149 + DENDRO_85;
const double DENDRO_151 = DENDRO_150 * DENDRO_98;
const double DENDRO_152 = DENDRO_54 * DENDRO_88;
const double DENDRO_153 = -0.25 * DENDRO_152;
const double DENDRO_154 = 0.5 * DENDRO_80;
const double DENDRO_155 = DENDRO_154 * DENDRO_69;
const double DENDRO_156 = -2 * DENDRO_103 * DENDRO_144 + DENDRO_155;
const double DENDRO_157 = 0.5 * DENDRO_69;
const double DENDRO_158 = DENDRO_157 * DENDRO_91;
const double DENDRO_159 = DENDRO_130 * (DENDRO_103 * DENDRO_135 + DENDRO_158);
const double DENDRO_160 = 0.5 * DENDRO_65;
const double DENDRO_161 = DENDRO_160 * DENDRO_88;
const double DENDRO_162 = DENDRO_160 * DENDRO_98;
const double DENDRO_163 = grad2_0_0_chi[pp];
const double DENDRO_164 = -DENDRO_163;
const double DENDRO_165 = -DENDRO_25;
const double DENDRO_166 = -DENDRO_29;
const double DENDRO_167 =
    DENDRO_10 * DENDRO_69 + DENDRO_165 * DENDRO_61 + DENDRO_166 * DENDRO_65;
const double DENDRO_168 = grad_2_chi[pp];
const double DENDRO_169 = DENDRO_16 * DENDRO_168;
const double DENDRO_170 = -DENDRO_33;
const double DENDRO_171 = DENDRO_102 + DENDRO_170 * DENDRO_69;
const double DENDRO_172 = grad_1_chi[pp];
const double DENDRO_173 = DENDRO_16 * DENDRO_172;
const double DENDRO_174 = -DENDRO_37;
const double DENDRO_175 =
    DENDRO_165 * DENDRO_65 + DENDRO_174 * DENDRO_61 + DENDRO_20 * DENDRO_69;
const double DENDRO_176 = DENDRO_16 * DENDRO_6;
const double DENDRO_177 = 1.0 / chi[pp];
const double DENDRO_178 = 2 * DENDRO_177;
const double DENDRO_179 = grad_2_gt3[pp];
const double DENDRO_180 = grad_1_gt5[pp];
const double DENDRO_181 = DENDRO_10 * DENDRO_180 + DENDRO_135 * DENDRO_20;
const double DENDRO_182 = -DENDRO_179 * DENDRO_33 + DENDRO_181;
const double DENDRO_183 = DENDRO_25 * DENDRO_91;
const double DENDRO_184 = grad_2_gt5[pp];
const double DENDRO_185 = 0.5 * DENDRO_184;
const double DENDRO_186 = DENDRO_10 * DENDRO_185;
const double DENDRO_187 = 0.5 * DENDRO_180;
const double DENDRO_188 = grad_2_gt4[pp];
const double DENDRO_189 = 1.0 * DENDRO_188;
const double DENDRO_190 = -DENDRO_189;
const double DENDRO_191 = DENDRO_187 + DENDRO_190;
const double DENDRO_192 =
    -DENDRO_150 * DENDRO_20 + DENDRO_186 + DENDRO_191 * DENDRO_33;
const double DENDRO_193 = DENDRO_192 * DENDRO_29;
const double DENDRO_194 = grad_1_gt3[pp];
const double DENDRO_195 = 0.5 * DENDRO_194;
const double DENDRO_196 = grad_1_gt4[pp];
const double DENDRO_197 = 1.0 * DENDRO_196;
const double DENDRO_198 = 0.5 * DENDRO_179;
const double DENDRO_199 = DENDRO_197 - DENDRO_198;
const double DENDRO_200 =
    -DENDRO_10 * DENDRO_199 + DENDRO_144 * DENDRO_20 + DENDRO_195 * DENDRO_33;
const double DENDRO_201 = -DENDRO_200;
const double DENDRO_202 = DENDRO_201 * DENDRO_33;
const double DENDRO_203 = DENDRO_104 + DENDRO_193 + DENDRO_202;
const double DENDRO_204 = DENDRO_10 * DENDRO_182 - 1.0 * DENDRO_183 +
                          DENDRO_20 * DENDRO_80 - DENDRO_203;
const double DENDRO_205 = 2.0 * DENDRO_49;
const double DENDRO_206 = DENDRO_204 * DENDRO_205;
const double DENDRO_207 =
    -DENDRO_10 * DENDRO_179 + DENDRO_135 * DENDRO_25 + DENDRO_180 * DENDRO_29;
const double DENDRO_208 = -DENDRO_207;
const double DENDRO_209 = DENDRO_25 * DENDRO_88;
const double DENDRO_210 =
    DENDRO_10 * DENDRO_191 - DENDRO_150 * DENDRO_25 + DENDRO_185 * DENDRO_29;
const double DENDRO_211 = -DENDRO_210;
const double DENDRO_212 = DENDRO_211 * DENDRO_29;
const double DENDRO_213 = DENDRO_10 * DENDRO_195;
const double DENDRO_214 =
    DENDRO_144 * DENDRO_25 - DENDRO_199 * DENDRO_29 + DENDRO_213;
const double DENDRO_215 = DENDRO_214 * DENDRO_33;
const double DENDRO_216 = DENDRO_110 + DENDRO_212 + DENDRO_215;
const double DENDRO_217 = DENDRO_10 * DENDRO_208 + DENDRO_20 * DENDRO_98 -
                          1.0 * DENDRO_209 - DENDRO_216;
const double DENDRO_218 = DENDRO_205 * DENDRO_217;
const double DENDRO_219 =
    DENDRO_135 * DENDRO_37 - DENDRO_179 * DENDRO_20 + DENDRO_180 * DENDRO_25;
const double DENDRO_220 = -DENDRO_219;
const double DENDRO_221 = DENDRO_25 * DENDRO_47;
const double DENDRO_222 =
    -DENDRO_150 * DENDRO_37 + DENDRO_185 * DENDRO_25 + DENDRO_191 * DENDRO_20;
const double DENDRO_223 = -DENDRO_222;
const double DENDRO_224 = DENDRO_223 * DENDRO_29;
const double DENDRO_225 = DENDRO_195 * DENDRO_20;
const double DENDRO_226 =
    DENDRO_144 * DENDRO_37 - DENDRO_199 * DENDRO_25 + DENDRO_225;
const double DENDRO_227 = DENDRO_226 * DENDRO_33;
const double DENDRO_228 = DENDRO_224 + DENDRO_227 + DENDRO_72;
const double DENDRO_229 = DENDRO_10 * DENDRO_220 + DENDRO_20 * DENDRO_56 -
                          1.0 * DENDRO_221 - DENDRO_228;
const double DENDRO_230 = DENDRO_205 * DENDRO_229;
const double DENDRO_231 = grad2_2_2_chi[pp];
const double DENDRO_232 = (DENDRO_168 * DENDRO_168);
const double DENDRO_233 = 3 * DENDRO_177;
const double DENDRO_234 = grad2_1_1_chi[pp];
const double DENDRO_235 = (DENDRO_172 * DENDRO_172);
const double DENDRO_236 = grad2_1_2_chi[pp];
const double DENDRO_237 = DENDRO_168 * DENDRO_233;
const double DENDRO_238 = grad2_0_2_chi[pp];
const double DENDRO_239 = 2 * DENDRO_25;
const double DENDRO_240 = grad2_0_1_chi[pp];
const double DENDRO_241 = DENDRO_172 * DENDRO_6;
const double DENDRO_242 =
    -2 * DENDRO_10 * (-DENDRO_172 * DENDRO_237 + 2 * DENDRO_236) +
    2 * DENDRO_169 * DENDRO_217 + 2 * DENDRO_173 * DENDRO_204 +
    2 * DENDRO_176 * DENDRO_229 -
    2 * DENDRO_20 * (-DENDRO_233 * DENDRO_241 + 2 * DENDRO_240) +
    DENDRO_239 * (-DENDRO_237 * DENDRO_6 + 2 * DENDRO_238) +
    DENDRO_29 * (2 * DENDRO_231 - DENDRO_232 * DENDRO_233) +
    DENDRO_33 * (-DENDRO_233 * DENDRO_235 + 2 * DENDRO_234) +
    DENDRO_37 * (2 * DENDRO_163 - DENDRO_233 * DENDRO_7);
const double DENDRO_243 = DENDRO_16 * DENDRO_177;
const double DENDRO_244 = -DENDRO_242 * DENDRO_243;
const double DENDRO_245 =
    4 * DENDRO_0 * gt0[pp] + DENDRO_1 * DENDRO_2 -
    DENDRO_101 * DENDRO_81 * DENDRO_98 - DENDRO_105 * DENDRO_107 * DENDRO_110 -
    DENDRO_106 + DENDRO_114 * DENDRO_137 +
    DENDRO_114 * (DENDRO_112 + DENDRO_47 * DENDRO_68) +
    DENDRO_114 * (DENDRO_116 + DENDRO_56 * DENDRO_64) +
    DENDRO_114 * (-1.0 * DENDRO_145 - DENDRO_148) +
    DENDRO_114 * (-1.0 * DENDRO_151 - DENDRO_153) +
    DENDRO_114 * (DENDRO_136 * DENDRO_88 + DENDRO_138) +
    DENDRO_119 * DENDRO_121 + DENDRO_121 * (DENDRO_127 + DENDRO_128) +
    DENDRO_121 * (DENDRO_142 + DENDRO_40 * DENDRO_98) -
    DENDRO_124 * (DENDRO_122 + DENDRO_40 * DENDRO_88) -
    DENDRO_124 * (DENDRO_125 + DENDRO_126) -
    DENDRO_130 * (1.0 * DENDRO_126 + DENDRO_129) -
    DENDRO_130 * (-2 * DENDRO_109 * DENDRO_150 + DENDRO_161) +
    DENDRO_132 * DENDRO_156 + DENDRO_132 * (1.0 * DENDRO_128 + DENDRO_131) +
    DENDRO_132 * (DENDRO_109 * DENDRO_135 + DENDRO_162) - DENDRO_141 -
    DENDRO_159 -
    DENDRO_178 * (DENDRO_164 + DENDRO_167 * DENDRO_169 +
                  DENDRO_171 * DENDRO_173 + DENDRO_175 * DENDRO_176) -
    DENDRO_19 + DENDRO_206 * DENDRO_52 + DENDRO_218 * DENDRO_40 - DENDRO_23 +
    DENDRO_230 * DENDRO_60 + DENDRO_24 * DENDRO_27 + DENDRO_244 * gt0[pp] +
    DENDRO_28 * DENDRO_31 + DENDRO_3 * DENDRO_4 + DENDRO_32 * DENDRO_35 +
    DENDRO_36 * DENDRO_39 - DENDRO_48 * DENDRO_51 - DENDRO_57 * DENDRO_59 -
    DENDRO_60 * DENDRO_72 * DENDRO_73 + DENDRO_78 * DENDRO_80 * DENDRO_81 -
    DENDRO_8 + DENDRO_86 * DENDRO_88 * DENDRO_89 - DENDRO_96;
const double DENDRO_246 = (x * x);
const double DENDRO_247 = (z * z);
const double DENDRO_248 = DENDRO_246 * gt0[pp];
const double DENDRO_249 = (y * y);
const double DENDRO_250 = DENDRO_249 * gt3[pp];
const double DENDRO_251 = x * y;
const double DENDRO_252 = 2 * DENDRO_251 * gt1[pp];
const double DENDRO_253 = DENDRO_246 + DENDRO_249;
const double DENDRO_254 = gt2[pp] * x;
const double DENDRO_255 = DENDRO_253 * DENDRO_254;
const double DENDRO_256 = 2 * z;
const double DENDRO_257 = gt4[pp] * y;
const double DENDRO_258 = DENDRO_253 * DENDRO_257;
const double DENDRO_259 = DENDRO_247 * gt5[pp] + DENDRO_248 + DENDRO_250 +
                          DENDRO_252 + DENDRO_254 * DENDRO_256 +
                          DENDRO_256 * DENDRO_257;
const double DENDRO_260 = 1.0 / DENDRO_259;
const double DENDRO_261 = DENDRO_247 * DENDRO_254 + DENDRO_247 * DENDRO_257 +
                          DENDRO_248 * z + DENDRO_250 * z + DENDRO_252 * z -
                          DENDRO_253 * gt5[pp] * z - DENDRO_255 - DENDRO_258;
const double DENDRO_262 = DENDRO_247 * DENDRO_248 + DENDRO_247 * DENDRO_250 +
                          DENDRO_247 * DENDRO_252 +
                          (DENDRO_253 * DENDRO_253) * gt5[pp] -
                          DENDRO_255 * DENDRO_256 - DENDRO_256 * DENDRO_258 -
                          DENDRO_260 * (DENDRO_261 * DENDRO_261);
const double DENDRO_263 = 1.0 / DENDRO_262;
const double DENDRO_264 = DENDRO_260 * DENDRO_261;
const double DENDRO_265 = -DENDRO_264 + z;
const double DENDRO_266 = DENDRO_263 * (DENDRO_265 * DENDRO_265);
const double DENDRO_267 = DENDRO_251 * gt0[pp];
const double DENDRO_268 = DENDRO_249 * gt1[pp];
const double DENDRO_269 = -DENDRO_246 * gt1[pp] + DENDRO_267 + DENDRO_268 +
                          gt2[pp] * y * z - gt3[pp] * x * y - gt4[pp] * x * z;
const double DENDRO_270 = -DENDRO_269;
const double DENDRO_271 = DENDRO_253 + DENDRO_264 * z;
const double DENDRO_272 = -DENDRO_246 * DENDRO_265 * gt1[pp] +
                          DENDRO_265 * DENDRO_267 + DENDRO_265 * DENDRO_268 -
                          DENDRO_265 * gt3[pp] * x * y -
                          DENDRO_271 * gt2[pp] * y + DENDRO_271 * gt4[pp] * x;
const double DENDRO_273 = -DENDRO_272;
const double DENDRO_274 = -DENDRO_246 * gt3[pp] - DENDRO_249 * gt0[pp] +
                          DENDRO_252 + DENDRO_260 * (DENDRO_270 * DENDRO_270) +
                          DENDRO_263 * (DENDRO_273 * DENDRO_273);
const double DENDRO_275 = 1.0 / DENDRO_274;
const double DENDRO_276 = DENDRO_260 * DENDRO_270;
const double DENDRO_277 = DENDRO_265 * x;
const double DENDRO_278 = DENDRO_263 * DENDRO_273;
const double DENDRO_279 = DENDRO_276 * x + DENDRO_277 * DENDRO_278 + y;
const double DENDRO_280 = chi[pp] * (DENDRO_246 * DENDRO_266 +
                                     DENDRO_275 * (DENDRO_279 * DENDRO_279));
const double DENDRO_281 = -DENDRO_234;
const double DENDRO_282 = -DENDRO_144;
const double DENDRO_283 =
    DENDRO_165 * DENDRO_282 + DENDRO_166 * DENDRO_199 + DENDRO_213;
const double DENDRO_284 =
    DENDRO_10 * DENDRO_199 + DENDRO_170 * DENDRO_195 + DENDRO_20 * DENDRO_282;
const double DENDRO_285 =
    DENDRO_165 * DENDRO_199 + DENDRO_174 * DENDRO_282 + DENDRO_225;
const double DENDRO_286 = DENDRO_220 * DENDRO_69;
const double DENDRO_287 = DENDRO_135 * DENDRO_56;
const double DENDRO_288 = 0.25 * DENDRO_287;
const double DENDRO_289 = DENDRO_191 * DENDRO_98;
const double DENDRO_290 = DENDRO_208 * DENDRO_54;
const double DENDRO_291 = -0.25 * DENDRO_290;
const double DENDRO_292 = DENDRO_180 * DENDRO_98;
const double DENDRO_293 = 0.25 * DENDRO_292;
const double DENDRO_294 = 0.5 * DENDRO_45;
const double DENDRO_295 = DENDRO_220 * DENDRO_52;
const double DENDRO_296 = 0.25 * DENDRO_295;
const double DENDRO_297 = 0.5 * DENDRO_56;
const double DENDRO_298 = 0.5 * DENDRO_144;
const double DENDRO_299 = DENDRO_220 * DENDRO_298;
const double DENDRO_300 = 0.5 * DENDRO_199;
const double DENDRO_301 = DENDRO_208 * DENDRO_300;
const double DENDRO_302 = 2 * DENDRO_191 * DENDRO_214;
const double DENDRO_303 = DENDRO_182 * DENDRO_194;
const double DENDRO_304 = 0.25 * DENDRO_303;
const double DENDRO_305 = DENDRO_179 * DENDRO_201;
const double DENDRO_306 = DENDRO_300 * DENDRO_98;
const double DENDRO_307 = -DENDRO_144 * DENDRO_297;
const double DENDRO_308 = 2 * DENDRO_226 * DENDRO_69;
const double DENDRO_309 = DENDRO_194 * DENDRO_80;
const double DENDRO_310 = 0.25 * DENDRO_309;
const double DENDRO_311 = DENDRO_201 * DENDRO_53;
const double DENDRO_312 = DENDRO_135 * DENDRO_226;
const double DENDRO_313 = 2.0 * DENDRO_113;
const double DENDRO_314 = DENDRO_180 * DENDRO_214;
const double DENDRO_315 = DENDRO_214 * DENDRO_54;
const double DENDRO_316 = DENDRO_226 * DENDRO_52;
const double DENDRO_317 = 0.25 * DENDRO_180;
const double DENDRO_318 = DENDRO_190 + DENDRO_317;
const double DENDRO_319 = DENDRO_208 * DENDRO_89;
const double DENDRO_320 = -DENDRO_94;
const double DENDRO_321 = DENDRO_89 * (DENDRO_320 + DENDRO_92 + DENDRO_95);
const double DENDRO_322 = DENDRO_143 + DENDRO_76;
const double DENDRO_323 = DENDRO_105 * DENDRO_227;
const double DENDRO_324 = DENDRO_37 * DENDRO_49;
const double DENDRO_325 = 4 * DENDRO_324;
const double DENDRO_326 = DENDRO_325 * (-DENDRO_115 + DENDRO_67);
const double DENDRO_327 = 0.75 * DENDRO_42;
const double DENDRO_328 = DENDRO_325 * (DENDRO_100 + DENDRO_320 + DENDRO_327);
const double DENDRO_329 = grad_1_Gt0[pp];
const double DENDRO_330 = grad_1_Gt1[pp];
const double DENDRO_331 = grad_1_Gt2[pp];
const double DENDRO_332 = 4 * gt4[pp];
const double DENDRO_333 = 0.25 * DENDRO_179;
const double DENDRO_334 = DENDRO_333 * DENDRO_80;
const double DENDRO_335 = DENDRO_182 * DENDRO_74;
const double DENDRO_336 = DENDRO_179 * DENDRO_182;
const double DENDRO_337 = DENDRO_53 * DENDRO_80;
const double DENDRO_338 = 3.0 * DENDRO_324;
const double DENDRO_339 =
    -DENDRO_105 * DENDRO_215 * (DENDRO_197 + DENDRO_198) -
    DENDRO_130 * (DENDRO_143 * DENDRO_182 + DENDRO_334) -
    DENDRO_130 * (DENDRO_198 * DENDRO_80 + DENDRO_335) -
    DENDRO_18 * grad2_1_2_gt3[pp] + DENDRO_2 * DENDRO_329 -
    DENDRO_22 * grad2_0_1_gt3[pp] - DENDRO_235 * DENDRO_5 +
    DENDRO_27 * grad2_0_2_gt3[pp] + DENDRO_31 * grad2_2_2_gt3[pp] +
    4 * DENDRO_330 * gt3[pp] + DENDRO_331 * DENDRO_332 -
    DENDRO_336 * DENDRO_51 - DENDRO_337 * DENDRO_338 +
    DENDRO_35 * grad2_1_1_gt3[pp] + DENDRO_39 * grad2_0_0_gt3[pp];
const double DENDRO_340 =
    DENDRO_114 * (DENDRO_301 - DENDRO_302) +
    DENDRO_114 * (DENDRO_304 + 1.0 * DENDRO_305) +
    DENDRO_114 * (DENDRO_226 * DENDRO_45 - DENDRO_299) +
    DENDRO_121 * (DENDRO_309 + DENDRO_311) +
    DENDRO_121 * (DENDRO_316 + DENDRO_53 * DENDRO_56) +
    DENDRO_121 * (DENDRO_179 * DENDRO_98 + DENDRO_315) -
    DENDRO_130 * (DENDRO_286 + DENDRO_288) -
    DENDRO_130 * (-1.0 * DENDRO_289 - DENDRO_291) -
    DENDRO_130 * (DENDRO_296 + DENDRO_297 * DENDRO_45) -
    DENDRO_130 * (DENDRO_208 * DENDRO_294 + DENDRO_293) +
    DENDRO_132 * (DENDRO_307 + DENDRO_308) +
    DENDRO_132 * (DENDRO_310 + 1.0 * DENDRO_311) +
    DENDRO_132 * (DENDRO_214 * DENDRO_45 + DENDRO_306) -
    DENDRO_178 * (DENDRO_169 * DENDRO_283 + DENDRO_173 * DENDRO_284 +
                  DENDRO_176 * DENDRO_285 + DENDRO_281) +
    DENDRO_179 * DENDRO_218 - DENDRO_194 * DENDRO_202 * DENDRO_73 +
    DENDRO_194 * DENDRO_206 - DENDRO_220 * DENDRO_321 + DENDRO_230 * DENDRO_53 +
    DENDRO_244 * gt3[pp] + DENDRO_313 * (DENDRO_303 + DENDRO_305) +
    DENDRO_313 * (DENDRO_179 * DENDRO_208 + DENDRO_314) +
    DENDRO_313 * (DENDRO_220 * DENDRO_53 + DENDRO_312) +
    DENDRO_318 * DENDRO_319 - DENDRO_322 * DENDRO_323 - DENDRO_326 * DENDRO_56 -
    DENDRO_328 * DENDRO_98 + DENDRO_339;
const double DENDRO_341 = DENDRO_265 * y;
const double DENDRO_342 = DENDRO_276 * y + DENDRO_278 * DENDRO_341 - x;
const double DENDRO_343 = chi[pp] * (DENDRO_249 * DENDRO_266 +
                                     DENDRO_275 * (DENDRO_342 * DENDRO_342));
const double DENDRO_344 = DENDRO_263 * chi[pp];
const double DENDRO_345 = -DENDRO_274;
const double DENDRO_346 = 1.0 / DENDRO_345;
const double DENDRO_347 = chi[pp] * z;
const double DENDRO_348 = DENDRO_271 * DENDRO_344;
const double DENDRO_349 = -DENDRO_273 * DENDRO_348 + DENDRO_276 * DENDRO_347;
const double DENDRO_350 = -DENDRO_349;
const double DENDRO_351 = DENDRO_177 * DENDRO_346 * (DENDRO_350 * DENDRO_350) -
                          (DENDRO_271 * DENDRO_271) * DENDRO_344;
const double DENDRO_352 = -DENDRO_231;
const double DENDRO_353 = -DENDRO_191;
const double DENDRO_354 = -DENDRO_150;
const double DENDRO_355 =
    DENDRO_10 * DENDRO_353 + DENDRO_165 * DENDRO_354 + DENDRO_166 * DENDRO_185;
const double DENDRO_356 =
    DENDRO_170 * DENDRO_353 + DENDRO_186 + DENDRO_20 * DENDRO_354;
const double DENDRO_357 =
    DENDRO_165 * DENDRO_185 + DENDRO_174 * DENDRO_354 + DENDRO_20 * DENDRO_353;
const double DENDRO_358 = 0.5 * DENDRO_191;
const double DENDRO_359 = DENDRO_358 * DENDRO_91;
const double DENDRO_360 = 0.5 * DENDRO_47;
const double DENDRO_361 = -DENDRO_150 * DENDRO_360;
const double DENDRO_362 = 2 * DENDRO_65;
const double DENDRO_363 = DENDRO_184 * DENDRO_88;
const double DENDRO_364 = 0.25 * DENDRO_363;
const double DENDRO_365 = DENDRO_211 * DENDRO_41;
const double DENDRO_366 = 0.5 * DENDRO_150;
const double DENDRO_367 = DENDRO_220 * DENDRO_366;
const double DENDRO_368 = DENDRO_184 * DENDRO_208;
const double DENDRO_369 = 0.25 * DENDRO_368;
const double DENDRO_370 = DENDRO_180 * DENDRO_211;
const double DENDRO_371 = DENDRO_220 * DENDRO_65;
const double DENDRO_372 = DENDRO_135 * DENDRO_47;
const double DENDRO_373 = 0.25 * DENDRO_372;
const double DENDRO_374 = DENDRO_317 * DENDRO_88;
const double DENDRO_375 = DENDRO_208 * DENDRO_82;
const double DENDRO_376 = DENDRO_220 * DENDRO_40;
const double DENDRO_377 = 0.25 * DENDRO_376;
const double DENDRO_378 = DENDRO_220 * DENDRO_41;
const double DENDRO_379 = DENDRO_135 * DENDRO_223;
const double DENDRO_380 = DENDRO_41 * DENDRO_47;
const double DENDRO_381 = DENDRO_223 * DENDRO_40;
const double DENDRO_382 = DENDRO_149 + DENDRO_84;
const double DENDRO_383 = DENDRO_105 * DENDRO_224;
const double DENDRO_384 = DENDRO_187 + DENDRO_189;
const double DENDRO_385 = DENDRO_105 * DENDRO_193;
const double DENDRO_386 = -DENDRO_100;
const double DENDRO_387 = DENDRO_81 * (DENDRO_386 + DENDRO_92 + DENDRO_99);
const double DENDRO_388 = DENDRO_325 * (-DENDRO_111 + DENDRO_63);
const double DENDRO_389 = DENDRO_180 * DENDRO_208;
const double DENDRO_390 = DENDRO_41 * DENDRO_88;
const double DENDRO_391 = grad_2_Gt0[pp];
const double DENDRO_392 = grad_2_Gt1[pp];
const double DENDRO_393 = grad_2_Gt2[pp];
const double DENDRO_394 = -DENDRO_182 * DENDRO_358;
const double DENDRO_395 = DENDRO_199 * DENDRO_91;
const double DENDRO_396 = DENDRO_182 * DENDRO_45;
const double DENDRO_397 = 0.25 * DENDRO_396;
const double DENDRO_398 = DENDRO_179 * DENDRO_91;
const double DENDRO_399 = 0.25 * DENDRO_398;
const double DENDRO_400 = 0.5 * DENDRO_182;
const double DENDRO_401 = DENDRO_180 * DENDRO_182;
const double DENDRO_402 = DENDRO_179 * DENDRO_192;
const double DENDRO_403 = DENDRO_180 * DENDRO_91;
const double DENDRO_404 = DENDRO_192 * DENDRO_45;
const double DENDRO_405 =
    DENDRO_114 * (2 * DENDRO_192 * DENDRO_199 + DENDRO_394) -
    DENDRO_124 * (DENDRO_403 + DENDRO_404) +
    DENDRO_132 * (DENDRO_395 + DENDRO_397) +
    DENDRO_132 * (DENDRO_399 + DENDRO_400 * DENDRO_54) -
    DENDRO_18 * grad2_1_2_gt5[pp] -
    DENDRO_182 * DENDRO_81 * (DENDRO_197 - DENDRO_333) -
    DENDRO_22 * grad2_0_1_gt5[pp] - DENDRO_232 * DENDRO_5 +
    DENDRO_27 * grad2_0_2_gt5[pp] + DENDRO_31 * grad2_2_2_gt5[pp] +
    DENDRO_313 * (DENDRO_401 + DENDRO_402) -
    DENDRO_325 * DENDRO_91 * (DENDRO_327 + DENDRO_386 + DENDRO_94) +
    DENDRO_332 * DENDRO_392 + DENDRO_35 * grad2_1_1_gt5[pp] +
    DENDRO_39 * grad2_0_0_gt5[pp] + DENDRO_391 * DENDRO_4 +
    4 * DENDRO_393 * gt5[pp];
const double DENDRO_406 =
    DENDRO_114 * (DENDRO_369 + 1.0 * DENDRO_370) +
    DENDRO_114 * (DENDRO_223 * DENDRO_54 - DENDRO_367) -
    DENDRO_124 * (DENDRO_363 + DENDRO_365) -
    DENDRO_124 * (DENDRO_380 + DENDRO_381) -
    DENDRO_130 * (DENDRO_364 + 1.0 * DENDRO_365) -
    DENDRO_130 * (DENDRO_192 * DENDRO_54 - DENDRO_359) -
    DENDRO_130 * (DENDRO_223 * DENDRO_362 + DENDRO_361) +
    DENDRO_132 * (DENDRO_371 + DENDRO_373) +
    DENDRO_132 * (DENDRO_149 * DENDRO_208 + DENDRO_374) +
    DENDRO_132 * (DENDRO_187 * DENDRO_88 + DENDRO_375) +
    DENDRO_132 * (DENDRO_360 * DENDRO_54 + DENDRO_377) -
    DENDRO_178 * (DENDRO_169 * DENDRO_355 + DENDRO_173 * DENDRO_356 +
                  DENDRO_176 * DENDRO_357 + DENDRO_352) +
    DENDRO_180 * DENDRO_206 - DENDRO_184 * DENDRO_212 * DENDRO_73 +
    DENDRO_184 * DENDRO_218 - DENDRO_220 * DENDRO_387 + DENDRO_230 * DENDRO_41 +
    DENDRO_244 * gt5[pp] + DENDRO_313 * (DENDRO_368 + DENDRO_370) +
    DENDRO_313 * (DENDRO_378 + DENDRO_379) - DENDRO_338 * DENDRO_390 -
    DENDRO_382 * DENDRO_383 - DENDRO_384 * DENDRO_385 - DENDRO_388 * DENDRO_47 -
    DENDRO_389 * DENDRO_59 + DENDRO_405;
const double DENDRO_407 = pow(DENDRO_177 * DENDRO_259, -1.0 / 2.0);
const double DENDRO_408 = At0[pp] * x;
const double DENDRO_409 = At1[pp] * y + At2[pp] * z + DENDRO_408;
const double DENDRO_410 = DENDRO_277 * DENDRO_409;
const double DENDRO_411 = At1[pp] * x + At3[pp] * y + At4[pp] * z;
const double DENDRO_412 = DENDRO_341 * DENDRO_411;
const double DENDRO_413 = At2[pp] * x + At4[pp] * y + At5[pp] * z;
const double DENDRO_414 = DENDRO_271 * DENDRO_413;
const double DENDRO_415 = DENDRO_410 + DENDRO_412 - DENDRO_414;
const double DENDRO_416 = 2 * DENDRO_407;
const double DENDRO_417 = -DENDRO_168 * DENDRO_271 + DENDRO_172 * DENDRO_341 +
                          DENDRO_277 * DENDRO_6 + DENDRO_410 * DENDRO_416 +
                          DENDRO_412 * DENDRO_416 - DENDRO_414 * DENDRO_416;
const double DENDRO_418 = DENDRO_279 * DENDRO_409;
const double DENDRO_419 = DENDRO_342 * DENDRO_411;
const double DENDRO_420 = DENDRO_177 * DENDRO_349;
const double DENDRO_421 = DENDRO_413 * DENDRO_420 + DENDRO_418 + DENDRO_419;
const double DENDRO_422 = DENDRO_168 * DENDRO_420 + DENDRO_172 * DENDRO_342 +
                          DENDRO_178 * DENDRO_349 * DENDRO_407 * DENDRO_413 +
                          DENDRO_279 * DENDRO_6 + DENDRO_416 * DENDRO_418 +
                          DENDRO_416 * DENDRO_419;
const double DENDRO_423 = DENDRO_346 * chi[pp];
const double DENDRO_424 =
    DENDRO_251 * DENDRO_266 + DENDRO_275 * DENDRO_279 * DENDRO_342;
const double DENDRO_425 = DENDRO_179 * DENDRO_88;
const double DENDRO_426 = DENDRO_135 * DENDRO_208;
const double DENDRO_427 = DENDRO_292 + DENDRO_426;
const double DENDRO_428 = DENDRO_47 * DENDRO_53;
const double DENDRO_429 = DENDRO_287 + DENDRO_295;
const double DENDRO_430 = DENDRO_208 * DENDRO_366;
const double DENDRO_431 = DENDRO_358 * DENDRO_88;
const double DENDRO_432 = DENDRO_430 + DENDRO_431;
const double DENDRO_433 = -0.25 * DENDRO_180 * DENDRO_88 + DENDRO_432;
const double DENDRO_434 = DENDRO_144 * DENDRO_201;
const double DENDRO_435 = -0.5 * DENDRO_194 * DENDRO_80 + DENDRO_434;
const double DENDRO_436 = DENDRO_136 * DENDRO_214 + DENDRO_333 * DENDRO_98;
const double DENDRO_437 = DENDRO_56 * DENDRO_74;
const double DENDRO_438 = DENDRO_226 * DENDRO_68;
const double DENDRO_439 = DENDRO_69 * DENDRO_71;
const double DENDRO_440 = DENDRO_109 * DENDRO_294 + DENDRO_111 * DENDRO_98;
const double DENDRO_441 = DENDRO_109 * DENDRO_191 + 0.5 * DENDRO_151;
const double DENDRO_442 = DENDRO_153 + DENDRO_441;
const double DENDRO_443 = DENDRO_111 * DENDRO_208;
const double DENDRO_444 = DENDRO_109 * DENDRO_187;
const double DENDRO_445 = DENDRO_45 * DENDRO_88;
const double DENDRO_446 = 0.25 * DENDRO_445;
const double DENDRO_447 = DENDRO_294 * DENDRO_71;
const double DENDRO_448 = DENDRO_112 + DENDRO_116;
const double DENDRO_449 = DENDRO_360 * DENDRO_69;
const double DENDRO_450 = 0.25 * DENDRO_220;
const double DENDRO_451 = DENDRO_450 * DENDRO_60;
const double DENDRO_452 = DENDRO_136 * DENDRO_71 + DENDRO_451;
const double DENDRO_453 = DENDRO_150 * DENDRO_214 + 0.5 * DENDRO_289;
const double DENDRO_454 = -0.5 * DENDRO_199 * DENDRO_88 + DENDRO_453;
const double DENDRO_455 = -DENDRO_144 * DENDRO_360;
const double DENDRO_456 = DENDRO_226 * DENDRO_64;
const double DENDRO_457 = 0.25 * DENDRO_45;
const double DENDRO_458 = DENDRO_456 + DENDRO_457 * DENDRO_56;
const double DENDRO_459 = DENDRO_136 * DENDRO_201;
const double DENDRO_460 = 0.25 * DENDRO_91;
const double DENDRO_461 = DENDRO_194 * DENDRO_460;
const double DENDRO_462 = DENDRO_334 + DENDRO_461;
const double DENDRO_463 = 0.25 * DENDRO_98;
const double DENDRO_464 = DENDRO_214 * DENDRO_64 + DENDRO_45 * DENDRO_463;
const double DENDRO_465 = DENDRO_226 * DENDRO_61;
const double DENDRO_466 = DENDRO_297 * DENDRO_69 + DENDRO_465;
const double DENDRO_467 = 0.25 * DENDRO_54;
const double DENDRO_468 = DENDRO_467 * DENDRO_98;
const double DENDRO_469 = DENDRO_109 * DENDRO_198 + DENDRO_135 * DENDRO_463;
const double DENDRO_470 = DENDRO_103 * DENDRO_195;
const double DENDRO_471 = DENDRO_470 + DENDRO_74 * DENDRO_80;
const double DENDRO_472 = DENDRO_45 * DENDRO_47;
const double DENDRO_473 = 1.0 * DENDRO_50;
const double DENDRO_474 = -DENDRO_240;
const double DENDRO_475 =
    DENDRO_10 * DENDRO_53 + DENDRO_165 * DENDRO_52 + DENDRO_166 * DENDRO_54;
const double DENDRO_476 = 0.5 * DENDRO_169;
const double DENDRO_477 = DENDRO_170 * DENDRO_53 + DENDRO_79;
const double DENDRO_478 = 0.5 * DENDRO_173;
const double DENDRO_479 =
    DENDRO_165 * DENDRO_54 + DENDRO_174 * DENDRO_52 + DENDRO_20 * DENDRO_53;
const double DENDRO_480 = 0.5 * DENDRO_176;
const double DENDRO_481 = 2.0 * DENDRO_329 * gt0[pp];
const double DENDRO_482 = 2.0 * gt1[pp];
const double DENDRO_483 = DENDRO_0 * DENDRO_482;
const double DENDRO_484 = DENDRO_330 * DENDRO_482;
const double DENDRO_485 = 2.0 * gt2[pp];
const double DENDRO_486 = DENDRO_331 * DENDRO_485;
const double DENDRO_487 = 2.0 * DENDRO_1 * gt3[pp];
const double DENDRO_488 = 2.0 * gt4[pp];
const double DENDRO_489 = DENDRO_3 * DENDRO_488;
const double DENDRO_490 = -DENDRO_241 * DENDRO_5;
const double DENDRO_491 = DENDRO_27 * grad2_0_2_gt1[pp];
const double DENDRO_492 = DENDRO_31 * grad2_2_2_gt1[pp];
const double DENDRO_493 = DENDRO_35 * grad2_1_1_gt1[pp];
const double DENDRO_494 = DENDRO_39 * grad2_0_0_gt1[pp];
const double DENDRO_495 = -DENDRO_18 * grad2_1_2_gt1[pp];
const double DENDRO_496 = -DENDRO_22 * grad2_0_1_gt1[pp];
const double DENDRO_497 =
    -DENDRO_178 * (DENDRO_474 + DENDRO_475 * DENDRO_476 +
                   DENDRO_477 * DENDRO_478 + DENDRO_479 * DENDRO_480) +
    DENDRO_206 * DENDRO_75 + DENDRO_218 * DENDRO_43 + DENDRO_230 * DENDRO_66 +
    DENDRO_244 * gt1[pp] + DENDRO_481 + DENDRO_483 + DENDRO_484 + DENDRO_486 +
    DENDRO_487 + DENDRO_489 + DENDRO_490 + DENDRO_491 + DENDRO_492 +
    DENDRO_493 + DENDRO_494 + DENDRO_495 + DENDRO_496;
const double DENDRO_498 = DENDRO_135 * DENDRO_182;
const double DENDRO_499 = 0.25 * DENDRO_498;
const double DENDRO_500 = 0.25 * DENDRO_117;
const double DENDRO_501 = DENDRO_103 * DENDRO_198;
const double DENDRO_502 = 0.25 * DENDRO_135;
const double DENDRO_503 = DENDRO_501 + DENDRO_502 * DENDRO_80;
const double DENDRO_504 = DENDRO_115 * DENDRO_182 + DENDRO_501;
const double DENDRO_505 = -DENDRO_144 * DENDRO_400;
const double DENDRO_506 = DENDRO_334 + DENDRO_505;
const double DENDRO_507 = -DENDRO_144 * DENDRO_154;
const double DENDRO_508 = DENDRO_114 * (DENDRO_461 + DENDRO_506) -
                          DENDRO_130 * (DENDRO_134 + DENDRO_503) -
                          DENDRO_130 * (DENDRO_134 + DENDRO_504) +
                          DENDRO_132 * (DENDRO_471 + DENDRO_507) -
                          DENDRO_325 * (1.0 * DENDRO_118 + DENDRO_500) -
                          DENDRO_89 * (0.5 * DENDRO_398 + DENDRO_499);
const double DENDRO_509 =
    DENDRO_113 * (DENDRO_425 + DENDRO_427) +
    DENDRO_113 * (DENDRO_428 + DENDRO_429) - DENDRO_114 * DENDRO_454 +
    DENDRO_114 * (DENDRO_455 + DENDRO_458) +
    DENDRO_114 * (DENDRO_459 + DENDRO_462) +
    DENDRO_121 * (DENDRO_53 * DENDRO_71 + DENDRO_57) + DENDRO_130 * DENDRO_442 -
    DENDRO_130 * (DENDRO_447 + DENDRO_448) -
    DENDRO_130 * (DENDRO_449 + DENDRO_452) -
    DENDRO_130 * (DENDRO_443 + DENDRO_444 + DENDRO_446) +
    DENDRO_132 * (DENDRO_468 + DENDRO_469) +
    DENDRO_132 * (DENDRO_109 * DENDRO_199 + DENDRO_464) +
    DENDRO_132 * (-DENDRO_144 * DENDRO_71 + DENDRO_466) +
    DENDRO_132 * (DENDRO_201 * DENDRO_68 + DENDRO_471) -
    DENDRO_325 * (0.5 * DENDRO_142 + DENDRO_440) -
    DENDRO_325 * (DENDRO_131 + DENDRO_439 + DENDRO_68 * DENDRO_71) +
    DENDRO_433 * DENDRO_89 + DENDRO_435 * DENDRO_81 -
    DENDRO_473 * (DENDRO_372 + DENDRO_376 + DENDRO_472) + DENDRO_497 +
    DENDRO_508 - DENDRO_81 * (DENDRO_306 + DENDRO_436) -
    DENDRO_81 * (DENDRO_307 + DENDRO_437 + DENDRO_438);
const double DENDRO_510 = 6 * chi[pp];
const double DENDRO_511 =
    DENDRO_275 * DENDRO_279 * DENDRO_349 - DENDRO_277 * DENDRO_348;
const double DENDRO_512 = DENDRO_180 * DENDRO_80;
const double DENDRO_513 = DENDRO_398 + DENDRO_512;
const double DENDRO_514 = DENDRO_498 + DENDRO_513;
const double DENDRO_515 = DENDRO_41 * DENDRO_56;
const double DENDRO_516 = DENDRO_376 + DENDRO_515;
const double DENDRO_517 = DENDRO_150 * DENDRO_211;
const double DENDRO_518 = -0.5 * DENDRO_184 * DENDRO_88 + DENDRO_517;
const double DENDRO_519 = DENDRO_136 * DENDRO_192 + 0.25 * DENDRO_403;
const double DENDRO_520 = DENDRO_89 * (-DENDRO_359 + DENDRO_519);
const double DENDRO_521 = 0.25 * DENDRO_380;
const double DENDRO_522 = 0.25 * DENDRO_426;
const double DENDRO_523 = DENDRO_300 * DENDRO_80;
const double DENDRO_524 = DENDRO_81 * (DENDRO_506 + DENDRO_523);
const double DENDRO_525 = DENDRO_111 * DENDRO_88;
const double DENDRO_526 = DENDRO_65 * DENDRO_71;
const double DENDRO_527 = 0.5 * DENDRO_54;
const double DENDRO_528 = DENDRO_103 * DENDRO_527 + 0.25 * DENDRO_139;
const double DENDRO_529 = DENDRO_325 * (0.5 * DENDRO_140 + DENDRO_528);
const double DENDRO_530 = DENDRO_223 * DENDRO_61;
const double DENDRO_531 = DENDRO_360 * DENDRO_65 + DENDRO_530;
const double DENDRO_532 = DENDRO_192 * DENDRO_68 + DENDRO_460 * DENDRO_54;
const double DENDRO_533 = DENDRO_130 * (-DENDRO_103 * DENDRO_191 + DENDRO_532);
const double DENDRO_534 = DENDRO_45 * DENDRO_460;
const double DENDRO_535 = DENDRO_103 * DENDRO_187 + DENDRO_135 * DENDRO_460;
const double DENDRO_536 = DENDRO_130 * (DENDRO_534 + DENDRO_535);
const double DENDRO_537 = DENDRO_109 * DENDRO_185;
const double DENDRO_538 = DENDRO_537 + DENDRO_82 * DENDRO_88;
const double DENDRO_539 = -DENDRO_366 * DENDRO_88;
const double DENDRO_540 = DENDRO_358 * DENDRO_80;
const double DENDRO_541 = DENDRO_144 * DENDRO_192;
const double DENDRO_542 =
    -0.5 * DENDRO_199 * DENDRO_91 + DENDRO_540 + DENDRO_541;
const double DENDRO_543 = -DENDRO_366 * DENDRO_56;
const double DENDRO_544 = DENDRO_223 * DENDRO_68;
const double DENDRO_545 = DENDRO_467 * DENDRO_47 + DENDRO_544;
const double DENDRO_546 = DENDRO_136 * DENDRO_211;
const double DENDRO_547 = DENDRO_184 * DENDRO_463;
const double DENDRO_548 = DENDRO_374 + DENDRO_547;
const double DENDRO_549 = -DENDRO_430;
const double DENDRO_550 = DENDRO_103 * DENDRO_199 - 0.5 * DENDRO_145;
const double DENDRO_551 = DENDRO_147 + DENDRO_550;
const double DENDRO_552 = DENDRO_502 * DENDRO_88;
const double DENDRO_553 = DENDRO_138 + DENDRO_444;
const double DENDRO_554 = DENDRO_54 * DENDRO_80;
const double DENDRO_555 = 0.25 * DENDRO_554;
const double DENDRO_556 = DENDRO_504 + DENDRO_555;
const double DENDRO_557 = DENDRO_527 * DENDRO_71;
const double DENDRO_558 = DENDRO_297 * DENDRO_65;
const double DENDRO_559 = DENDRO_54 * DENDRO_56;
const double DENDRO_560 = 1.0 * DENDRO_58;
const double DENDRO_561 = -DENDRO_238;
const double DENDRO_562 =
    DENDRO_10 * DENDRO_45 + DENDRO_165 * DENDRO_40 + DENDRO_166 * DENDRO_41;
const double DENDRO_563 = DENDRO_170 * DENDRO_45 + DENDRO_90;
const double DENDRO_564 =
    DENDRO_165 * DENDRO_41 + DENDRO_174 * DENDRO_40 + DENDRO_20 * DENDRO_45;
const double DENDRO_565 = 2.0 * DENDRO_391;
const double DENDRO_566 = 2.0 * DENDRO_392;
const double DENDRO_567 = 2.0 * gt5[pp];
const double DENDRO_568 = DENDRO_168 * DENDRO_5;
const double DENDRO_569 = DENDRO_568 * DENDRO_6;
const double DENDRO_570 = grad2_0_2_gt2[pp];
const double DENDRO_571 = grad2_2_2_gt2[pp];
const double DENDRO_572 = grad2_1_1_gt2[pp];
const double DENDRO_573 = grad2_0_0_gt2[pp];
const double DENDRO_574 = DENDRO_18 * grad2_1_2_gt2[pp];
const double DENDRO_575 = DENDRO_22 * grad2_0_1_gt2[pp];
const double DENDRO_576 =
    DENDRO_0 * DENDRO_485 + DENDRO_1 * DENDRO_488 -
    DENDRO_178 * (DENDRO_476 * DENDRO_562 + DENDRO_478 * DENDRO_563 +
                  DENDRO_480 * DENDRO_564 + DENDRO_561) +
    DENDRO_206 * DENDRO_44 + DENDRO_218 * DENDRO_83 + DENDRO_230 * DENDRO_62 +
    DENDRO_244 * gt2[pp] + DENDRO_27 * DENDRO_570 + DENDRO_3 * DENDRO_567 +
    DENDRO_31 * DENDRO_571 + DENDRO_35 * DENDRO_572 + DENDRO_39 * DENDRO_573 +
    DENDRO_393 * DENDRO_485 + DENDRO_565 * gt0[pp] + DENDRO_566 * gt1[pp] -
    DENDRO_569 - DENDRO_574 - DENDRO_575;
const double DENDRO_577 =
    DENDRO_113 * DENDRO_514 + DENDRO_113 * (DENDRO_372 + DENDRO_516) -
    DENDRO_114 * DENDRO_542 + DENDRO_114 * (DENDRO_543 + DENDRO_545) +
    DENDRO_114 * (DENDRO_546 + DENDRO_548) +
    DENDRO_114 * (DENDRO_548 + DENDRO_549) -
    DENDRO_124 * (DENDRO_41 * DENDRO_71 + DENDRO_48) -
    DENDRO_130 * (DENDRO_538 + DENDRO_539) -
    DENDRO_130 * (-DENDRO_150 * DENDRO_71 + DENDRO_531) -
    DENDRO_130 * (DENDRO_211 * DENDRO_64 + DENDRO_538) +
    DENDRO_132 * DENDRO_551 + DENDRO_132 * DENDRO_556 +
    DENDRO_132 * (DENDRO_443 + DENDRO_553) +
    DENDRO_132 * (DENDRO_448 + DENDRO_557) +
    DENDRO_132 * (DENDRO_452 + DENDRO_558) +
    DENDRO_132 * (DENDRO_552 + DENDRO_553) -
    DENDRO_325 * (1.0 * DENDRO_122 + DENDRO_525) -
    DENDRO_325 * (DENDRO_129 + DENDRO_526 + DENDRO_64 * DENDRO_71) +
    DENDRO_518 * DENDRO_89 - DENDRO_520 - DENDRO_524 - DENDRO_529 - DENDRO_533 -
    DENDRO_536 - DENDRO_560 * (DENDRO_429 + DENDRO_559) + DENDRO_576 -
    DENDRO_81 * (DENDRO_187 * DENDRO_98 + DENDRO_522) -
    DENDRO_89 * (DENDRO_223 * DENDRO_64 + DENDRO_361 + DENDRO_521);
const double DENDRO_578 =
    DENDRO_275 * DENDRO_342 * DENDRO_349 - DENDRO_341 * DENDRO_348;
const double DENDRO_579 = DENDRO_191 * DENDRO_211;
const double DENDRO_580 = -0.5 * DENDRO_184 * DENDRO_208 + DENDRO_579;
const double DENDRO_581 = DENDRO_223 * DENDRO_294 + 0.25 * DENDRO_378;
const double DENDRO_582 = DENDRO_208 * DENDRO_333;
const double DENDRO_583 = 1.0 * DENDRO_314;
const double DENDRO_584 = DENDRO_199 * DENDRO_201;
const double DENDRO_585 = DENDRO_220 * DENDRO_74;
const double DENDRO_586 = DENDRO_226 * DENDRO_527;
const double DENDRO_587 = 0.5 * DENDRO_312 + DENDRO_586;
const double DENDRO_588 = DENDRO_449 + DENDRO_558;
const double DENDRO_589 = DENDRO_223 * DENDRO_69 + 0.5 * DENDRO_371;
const double DENDRO_590 = DENDRO_211 * DENDRO_294;
const double DENDRO_591 = DENDRO_375 + DENDRO_547;
const double DENDRO_592 = -DENDRO_431;
const double DENDRO_593 = DENDRO_150 * DENDRO_226;
const double DENDRO_594 = DENDRO_143 * DENDRO_223 + DENDRO_220 * DENDRO_467;
const double DENDRO_595 = DENDRO_192 * DENDRO_195;
const double DENDRO_596 = DENDRO_182 * DENDRO_300 + DENDRO_595;
const double DENDRO_597 = DENDRO_220 * DENDRO_502;
const double DENDRO_598 = DENDRO_149 * DENDRO_226 + DENDRO_45 * DENDRO_450;
const double DENDRO_599 = DENDRO_185 * DENDRO_214;
const double DENDRO_600 = DENDRO_208 * DENDRO_317 + DENDRO_599;
const double DENDRO_601 = -DENDRO_208 * DENDRO_358;
const double DENDRO_602 = DENDRO_226 * DENDRO_65 + 0.5 * DENDRO_286;
const double DENDRO_603 = DENDRO_208 * DENDRO_457;
const double DENDRO_604 = DENDRO_149 * DENDRO_214;
const double DENDRO_605 = DENDRO_293 + DENDRO_604;
const double DENDRO_606 = 0.25 * DENDRO_559;
const double DENDRO_607 = 0.25 * DENDRO_428 + DENDRO_456;
const double DENDRO_608 = DENDRO_333 * DENDRO_88;
const double DENDRO_609 = DENDRO_201 * DENDRO_527;
const double DENDRO_610 = DENDRO_334 + DENDRO_335;
const double DENDRO_611 = DENDRO_201 * DENDRO_294;
const double DENDRO_612 = DENDRO_461 + DENDRO_523;
const double DENDRO_613 = 1.0 * DENDRO_123;
const double DENDRO_614 = DENDRO_565 * gt1[pp];
const double DENDRO_615 = DENDRO_329 * DENDRO_485;
const double DENDRO_616 = DENDRO_566 * gt3[pp];
const double DENDRO_617 = DENDRO_330 * DENDRO_488;
const double DENDRO_618 = DENDRO_393 * DENDRO_488;
const double DENDRO_619 = DENDRO_331 * DENDRO_567;
const double DENDRO_620 = -DENDRO_172 * DENDRO_568;
const double DENDRO_621 = 0.25 * DENDRO_401;
const double DENDRO_622 = DENDRO_143 * DENDRO_192;
const double DENDRO_623 = DENDRO_182 * DENDRO_467 + DENDRO_622;
const double DENDRO_624 = DENDRO_27 * grad2_0_2_gt4[pp];
const double DENDRO_625 = DENDRO_31 * grad2_2_2_gt4[pp];
const double DENDRO_626 = DENDRO_35 * grad2_1_1_gt4[pp];
const double DENDRO_627 = DENDRO_39 * grad2_0_0_gt4[pp];
const double DENDRO_628 = DENDRO_133 + DENDRO_554;
const double DENDRO_629 = 1.0 * DENDRO_324;
const double DENDRO_630 = -DENDRO_18 * grad2_1_2_gt4[pp];
const double DENDRO_631 = -DENDRO_22 * grad2_0_1_gt4[pp];
const double DENDRO_632 =
    -DENDRO_130 * (-DENDRO_540 + DENDRO_623) -
    DENDRO_613 * (DENDRO_396 + DENDRO_513) + DENDRO_614 + DENDRO_615 +
    DENDRO_616 + DENDRO_617 + DENDRO_618 + DENDRO_619 + DENDRO_620 +
    DENDRO_624 + DENDRO_625 + DENDRO_626 + DENDRO_627 -
    DENDRO_629 * (DENDRO_146 + DENDRO_628) + DENDRO_630 + DENDRO_631 -
    DENDRO_89 * (DENDRO_192 * DENDRO_198 + DENDRO_394 + DENDRO_621);
const double DENDRO_633 = -DENDRO_236;
const double DENDRO_634 =
    DENDRO_10 * DENDRO_179 + DENDRO_135 * DENDRO_165 + DENDRO_166 * DENDRO_180;
const double DENDRO_635 = DENDRO_170 * DENDRO_179 + DENDRO_181;
const double DENDRO_636 =
    DENDRO_135 * DENDRO_174 + DENDRO_165 * DENDRO_180 + DENDRO_179 * DENDRO_20;
const double DENDRO_637 =
    -DENDRO_178 * (DENDRO_476 * DENDRO_634 + DENDRO_478 * DENDRO_635 +
                   DENDRO_480 * DENDRO_636 + DENDRO_633) +
    DENDRO_188 * DENDRO_218 + DENDRO_196 * DENDRO_206 + DENDRO_230 * DENDRO_42 +
    DENDRO_244 * gt4[pp];
const double DENDRO_638 =
    DENDRO_114 * (-DENDRO_593 + DENDRO_594) +
    DENDRO_114 * (DENDRO_597 + DENDRO_598) +
    DENDRO_114 * (DENDRO_600 + DENDRO_601) +
    DENDRO_114 * (-DENDRO_191 * DENDRO_201 + DENDRO_596) +
    DENDRO_114 * (DENDRO_198 * DENDRO_211 + DENDRO_600) -
    DENDRO_130 * (DENDRO_543 + DENDRO_589) -
    DENDRO_130 * (DENDRO_590 + DENDRO_591) -
    DENDRO_130 * (DENDRO_591 + DENDRO_592) +
    DENDRO_132 * (DENDRO_288 + DENDRO_602) +
    DENDRO_132 * (DENDRO_603 + DENDRO_605) +
    DENDRO_132 * (DENDRO_605 + DENDRO_608) +
    DENDRO_132 * (DENDRO_606 + DENDRO_607) +
    DENDRO_132 * (DENDRO_609 + DENDRO_610) +
    DENDRO_132 * (DENDRO_611 + DENDRO_612) +
    DENDRO_313 * (DENDRO_180 * DENDRO_201 + DENDRO_336) -
    DENDRO_325 * (DENDRO_112 + DENDRO_588) -
    DENDRO_325 * (DENDRO_149 * DENDRO_98 + DENDRO_446) +
    DENDRO_580 * DENDRO_89 - DENDRO_613 * (DENDRO_472 + DENDRO_516) +
    DENDRO_632 + DENDRO_637 - DENDRO_81 * (DENDRO_582 + DENDRO_583) -
    DENDRO_81 * (DENDRO_585 + DENDRO_587) -
    DENDRO_81 * (DENDRO_198 * DENDRO_201 + DENDRO_304 + DENDRO_584) -
    DENDRO_89 * (-DENDRO_367 + DENDRO_581);
const double DENDRO_639 = 2 * chi[pp];
const double DENDRO_640 = 2 * At2[pp];
const double DENDRO_641 = 2 * At4[pp];
const double DENDRO_642 = 3 * DENDRO_407;
const double DENDRO_643 = DENDRO_260 * DENDRO_510;
const double DENDRO_644 =
    2 * DENDRO_5 *
    (DENDRO_168 * DENDRO_642 * z + DENDRO_172 * DENDRO_642 * y +
     6 * DENDRO_260 * DENDRO_347 * DENDRO_413 + DENDRO_409 * DENDRO_643 * x +
     DENDRO_411 * DENDRO_643 * y + DENDRO_6 * DENDRO_642 * x +
     DENDRO_639 * K[pp]);
const double DENDRO_645 = grad_0_At2[pp];
const double DENDRO_646 = At1[pp] * DENDRO_16;
const double DENDRO_647 = 0.5 * DENDRO_91;
const double DENDRO_648 = At0[pp] * DENDRO_16;
const double DENDRO_649 = At2[pp] * DENDRO_16;
const double DENDRO_650 = 0.5 * DENDRO_649;
const double DENDRO_651 = DENDRO_360 * DENDRO_648 + DENDRO_645 +
                          DENDRO_646 * DENDRO_647 + DENDRO_650 * DENDRO_88;
const double DENDRO_652 = grad_0_At4[pp];
const double DENDRO_653 = 0.5 * DENDRO_220;
const double DENDRO_654 = DENDRO_208 * DENDRO_650 + DENDRO_400 * DENDRO_646 +
                          DENDRO_648 * DENDRO_653 + DENDRO_652;
const double DENDRO_655 = grad_0_At5[pp];
const double DENDRO_656 = DENDRO_192 * DENDRO_646 + DENDRO_211 * DENDRO_649 +
                          DENDRO_223 * DENDRO_648 + DENDRO_655;
const double DENDRO_657 = -DENDRO_271 * DENDRO_656 + DENDRO_277 * DENDRO_651 +
                          DENDRO_341 * DENDRO_654;
const double DENDRO_658 = DENDRO_265 * z + DENDRO_271;
const double DENDRO_659 = DENDRO_658 * x;
const double DENDRO_660 = grad_1_At2[pp];
const double DENDRO_661 = At3[pp] * DENDRO_16;
const double DENDRO_662 = At4[pp] * DENDRO_16;
const double DENDRO_663 = 0.5 * DENDRO_662;
const double DENDRO_664 = DENDRO_360 * DENDRO_646 + DENDRO_647 * DENDRO_661 +
                          DENDRO_660 + DENDRO_663 * DENDRO_88;
const double DENDRO_665 = grad_1_At4[pp];
const double DENDRO_666 = DENDRO_208 * DENDRO_663 + DENDRO_400 * DENDRO_661 +
                          DENDRO_646 * DENDRO_653 + DENDRO_665;
const double DENDRO_667 = grad_1_At5[pp];
const double DENDRO_668 = DENDRO_192 * DENDRO_661 + DENDRO_211 * DENDRO_662 +
                          DENDRO_223 * DENDRO_646 + DENDRO_667;
const double DENDRO_669 = -DENDRO_271 * DENDRO_668 + DENDRO_277 * DENDRO_664 +
                          DENDRO_341 * DENDRO_666;
const double DENDRO_670 = DENDRO_658 * y;
const double DENDRO_671 = grad_2_At1[pp];
const double DENDRO_672 = At5[pp] * DENDRO_16;
const double DENDRO_673 = 0.5 * DENDRO_672;
const double DENDRO_674 = DENDRO_154 * DENDRO_662 + DENDRO_297 * DENDRO_649 +
                          DENDRO_671 + DENDRO_673 * DENDRO_98;
const double DENDRO_675 = grad_2_At2[pp];
const double DENDRO_676 = DENDRO_360 * DENDRO_649 + DENDRO_647 * DENDRO_662 +
                          DENDRO_673 * DENDRO_88 + DENDRO_675;
const double DENDRO_677 = grad_2_At0[pp];
const double DENDRO_678 = DENDRO_103 * DENDRO_662 + DENDRO_109 * DENDRO_672 +
                          DENDRO_649 * DENDRO_71 + DENDRO_677;
const double DENDRO_679 = -DENDRO_271 * DENDRO_676 + DENDRO_277 * DENDRO_678 +
                          DENDRO_341 * DENDRO_674;
const double DENDRO_680 = grad_2_At4[pp];
const double DENDRO_681 = DENDRO_208 * DENDRO_673 + DENDRO_220 * DENDRO_650 +
                          DENDRO_400 * DENDRO_662 + DENDRO_680;
const double DENDRO_682 = grad_2_At3[pp];
const double DENDRO_683 = DENDRO_201 * DENDRO_662 + DENDRO_214 * DENDRO_672 +
                          DENDRO_226 * DENDRO_649 + DENDRO_682;
const double DENDRO_684 = -DENDRO_271 * DENDRO_681 + DENDRO_277 * DENDRO_674 +
                          DENDRO_341 * DENDRO_683;
const double DENDRO_685 = DENDRO_279 * y;
const double DENDRO_686 = -DENDRO_342;
const double DENDRO_687 = DENDRO_685 + DENDRO_686 * x;
const double DENDRO_688 = grad_1_At1[pp];
const double DENDRO_689 = DENDRO_154 * DENDRO_661 + DENDRO_297 * DENDRO_646 +
                          DENDRO_663 * DENDRO_98 + DENDRO_688;
const double DENDRO_690 = grad_1_At0[pp];
const double DENDRO_691 = DENDRO_103 * DENDRO_661 + DENDRO_109 * DENDRO_662 +
                          DENDRO_646 * DENDRO_71 + DENDRO_690;
const double DENDRO_692 = DENDRO_342 * x;
const double DENDRO_693 = DENDRO_279 * y - DENDRO_692;
const double DENDRO_694 = grad_0_At1[pp];
const double DENDRO_695 = DENDRO_154 * DENDRO_646 + DENDRO_297 * DENDRO_648 +
                          DENDRO_650 * DENDRO_98 + DENDRO_694;
const double DENDRO_696 = grad_0_At3[pp];
const double DENDRO_697 = DENDRO_201 * DENDRO_646 + DENDRO_214 * DENDRO_649 +
                          DENDRO_226 * DENDRO_648 + DENDRO_696;
const double DENDRO_698 = DENDRO_279 * z;
const double DENDRO_699 = DENDRO_177 * DENDRO_350 * x + DENDRO_698;
const double DENDRO_700 =
    DENDRO_279 * DENDRO_678 + DENDRO_342 * DENDRO_674 + DENDRO_420 * DENDRO_676;
const double DENDRO_701 = -DENDRO_420 * x + DENDRO_698;
const double DENDRO_702 =
    DENDRO_279 * DENDRO_651 + DENDRO_342 * DENDRO_654 + DENDRO_420 * DENDRO_656;
const double DENDRO_703 = DENDRO_177 * DENDRO_350 * y - DENDRO_686 * z;
const double DENDRO_704 =
    DENDRO_279 * DENDRO_674 + DENDRO_342 * DENDRO_683 + DENDRO_420 * DENDRO_681;
const double DENDRO_705 = DENDRO_342 * z - DENDRO_420 * y;
const double DENDRO_706 =
    DENDRO_279 * DENDRO_664 + DENDRO_342 * DENDRO_666 + DENDRO_420 * DENDRO_668;
const double DENDRO_707 = DENDRO_260 * DENDRO_269;
const double DENDRO_708 = DENDRO_263 * DENDRO_272;
const double DENDRO_709 = -DENDRO_277 * DENDRO_708 - DENDRO_707 * x + y;
const double DENDRO_710 = DENDRO_111 * DENDRO_55;
const double DENDRO_711 = DENDRO_115 * DENDRO_46;
const double DENDRO_712 = DENDRO_108 * DENDRO_41;
const double DENDRO_713 = DENDRO_46 * DENDRO_60;
const double DENDRO_714 = DENDRO_40 * DENDRO_70;
const double DENDRO_715 = DENDRO_55 * DENDRO_60;
const double DENDRO_716 = DENDRO_52 * DENDRO_70;
const double DENDRO_717 = 0.25 * DENDRO_713;
const double DENDRO_718 = 0.25 * DENDRO_715;
const double DENDRO_719 = DENDRO_82 * DENDRO_97;
const double DENDRO_720 = DENDRO_108 * DENDRO_54;
const double DENDRO_721 =
    DENDRO_205 * (-1.0 * DENDRO_10 * DENDRO_182 + DENDRO_183 -
                  1.0 * DENDRO_20 * DENDRO_80 + DENDRO_203);
const double DENDRO_722 =
    DENDRO_205 * (-1.0 * DENDRO_10 * DENDRO_208 - 1.0 * DENDRO_20 * DENDRO_98 +
                  DENDRO_209 + DENDRO_216);
const double DENDRO_723 =
    DENDRO_205 * (-1.0 * DENDRO_10 * DENDRO_220 - 1.0 * DENDRO_20 * DENDRO_56 +
                  DENDRO_221 + DENDRO_228);
const double DENDRO_724 = DENDRO_242 * DENDRO_243;
const double DENDRO_725 = DENDRO_341 * DENDRO_708 + DENDRO_707 * y + x;
const double DENDRO_726 = -0.25 * DENDRO_182 * DENDRO_194;
const double DENDRO_727 = DENDRO_179 * DENDRO_200;
const double DENDRO_728 = 1.0 * DENDRO_219;
const double DENDRO_729 = DENDRO_135 * DENDRO_55;
const double DENDRO_730 = 0.25 * DENDRO_729;
const double DENDRO_731 = DENDRO_317 * DENDRO_97;
const double DENDRO_732 = DENDRO_300 * DENDRO_97;
const double DENDRO_733 = DENDRO_298 * DENDRO_55;
const double DENDRO_734 = DENDRO_200 * DENDRO_53;
const double DENDRO_735 = DENDRO_184 * DENDRO_207;
const double DENDRO_736 = DENDRO_180 * DENDRO_210;
const double DENDRO_737 = -0.5 * DENDRO_150 * DENDRO_46;
const double DENDRO_738 = DENDRO_135 * DENDRO_46;
const double DENDRO_739 = DENDRO_317 * DENDRO_87;
const double DENDRO_740 = DENDRO_207 * DENDRO_82;
const double DENDRO_741 = DENDRO_184 * DENDRO_87;
const double DENDRO_742 = DENDRO_210 * DENDRO_41;
const double DENDRO_743 =
    DENDRO_263 * DENDRO_271 * DENDRO_272 * chi[pp] - DENDRO_347 * DENDRO_707;
const double DENDRO_744 = DENDRO_265 * (DENDRO_685 + DENDRO_692);
const double DENDRO_745 =
    DENDRO_177 * DENDRO_265 * DENDRO_349 * x - DENDRO_271 * DENDRO_279;
const double DENDRO_746 =
    DENDRO_177 * DENDRO_265 * DENDRO_349 * y - DENDRO_271 * DENDRO_342;
const double DENDRO_747 = DENDRO_219 * DENDRO_40;
const double DENDRO_748 = DENDRO_45 * DENDRO_46;
const double DENDRO_749 = DENDRO_111 * DENDRO_207;
const double DENDRO_750 = DENDRO_108 * DENDRO_187;
const double DENDRO_751 = DENDRO_457 * DENDRO_87;
const double DENDRO_752 = DENDRO_157 * DENDRO_46;
const double DENDRO_753 =
    DENDRO_136 * DENDRO_70 + 0.25 * DENDRO_219 * DENDRO_60;
const double DENDRO_754 = DENDRO_710 + DENDRO_711;
const double DENDRO_755 = DENDRO_135 * DENDRO_207;
const double DENDRO_756 = 1.0 * DENDRO_113;
const double DENDRO_757 = DENDRO_219 * DENDRO_52 + DENDRO_729;
const double DENDRO_758 = DENDRO_719 + DENDRO_750;
const double DENDRO_759 = 0.25 * DENDRO_184 * DENDRO_97;
const double DENDRO_760 = DENDRO_739 + DENDRO_759;
const double DENDRO_761 = DENDRO_41 * DENDRO_55 + DENDRO_747;
const double DENDRO_762 = DENDRO_108 * DENDRO_185 + DENDRO_82 * DENDRO_87;
const double DENDRO_763 = DENDRO_54 * DENDRO_55;
const double DENDRO_764 = -0.5 * DENDRO_150 * DENDRO_55;
const double DENDRO_765 = DENDRO_160 * DENDRO_55;
const double DENDRO_766 = DENDRO_740 + DENDRO_759;
const double DENDRO_767 = DENDRO_207 * DENDRO_317;
const double DENDRO_768 = -0.5 * DENDRO_214 * DENDRO_41 + DENDRO_731;
const double DENDRO_769 = (DENDRO_20 * DENDRO_20);
const double DENDRO_770 = (DENDRO_25 * DENDRO_25);
const double DENDRO_771 = At4[pp] * DENDRO_20;
const double DENDRO_772 = 2 * DENDRO_20;
const double DENDRO_773 = At1[pp] * DENDRO_37;
const double DENDRO_774 = At2[pp] * DENDRO_37;
const double DENDRO_775 = (DENDRO_10 * DENDRO_10);
const double DENDRO_776 = At2[pp] * DENDRO_10;
const double DENDRO_777 = 2 * DENDRO_10;
const double DENDRO_778 = At1[pp] * DENDRO_10;
const double DENDRO_779 = At2[pp] * DENDRO_29;
const double DENDRO_780 = At4[pp] * DENDRO_29;
const double DENDRO_781 = At0[pp] * DENDRO_37;
const double DENDRO_782 = DENDRO_20 * DENDRO_25;
const double DENDRO_783 = At5[pp] * DENDRO_25;
const double DENDRO_784 = At3[pp] * DENDRO_33;
const double DENDRO_785 = At5[pp] * DENDRO_10;
const double DENDRO_786 = (1.0 / 4.0) * chi[pp];
const double DENDRO_787 = DENDRO_208 * DENDRO_40;
const double DENDRO_788 = DENDRO_41 * DENDRO_98;
const double DENDRO_789 = DENDRO_152 + DENDRO_788;
const double DENDRO_790 = DENDRO_182 * DENDRO_52;
const double DENDRO_791 = 0.25 * DENDRO_48;
const double DENDRO_792 = 0.25 * DENDRO_512 + DENDRO_622;
const double DENDRO_793 = DENDRO_374 + DENDRO_375;
const double DENDRO_794 = DENDRO_211 * DENDRO_527 + DENDRO_547;
const double DENDRO_795 = 0.25 * DENDRO_515 + DENDRO_544;
const double DENDRO_796 = DENDRO_400 * DENDRO_69;
const double DENDRO_797 = DENDRO_116 + DENDRO_451;
const double DENDRO_798 = DENDRO_26 * DENDRO_786;
const double DENDRO_799 = 0.25 * DENDRO_472;
const double DENDRO_800 = 0.25 * DENDRO_336;
const double DENDRO_801 = DENDRO_112 + DENDRO_451;
const double DENDRO_802 = 0.25 * DENDRO_57;
const double DENDRO_803 = (3.0 / 2.0) * DENDRO_177;
const double DENDRO_804 = DENDRO_20 * DENDRO_803;
const double DENDRO_805 = DENDRO_173 * DENDRO_804;
const double DENDRO_806 = DENDRO_169 * DENDRO_803;
const double DENDRO_807 = DENDRO_176 * DENDRO_804;
const double DENDRO_808 = DENDRO_173 * DENDRO_803;
const double DENDRO_809 = -DENDRO_671;
const double DENDRO_810 = 0.5 * DENDRO_646;
const double DENDRO_811 = 0.5 * DENDRO_661;
const double DENDRO_812 =
    DENDRO_562 * DENDRO_663 + DENDRO_563 * DENDRO_811 + DENDRO_564 * DENDRO_810;
const double DENDRO_813 = -DENDRO_660;
const double DENDRO_814 =
    DENDRO_475 * DENDRO_673 + DENDRO_477 * DENDRO_663 + DENDRO_479 * DENDRO_650;
const double DENDRO_815 = 0.5 * DENDRO_648;
const double DENDRO_816 =
    DENDRO_475 * DENDRO_650 + DENDRO_477 * DENDRO_810 + DENDRO_479 * DENDRO_815;
const double DENDRO_817 =
    DENDRO_562 * DENDRO_650 + DENDRO_563 * DENDRO_810 + DENDRO_564 * DENDRO_815;
const double DENDRO_818 = -DENDRO_675;
const double DENDRO_819 =
    DENDRO_562 * DENDRO_673 + DENDRO_563 * DENDRO_663 + DENDRO_564 * DENDRO_650;
const double DENDRO_820 = -DENDRO_688;
const double DENDRO_821 =
    DENDRO_475 * DENDRO_663 + DENDRO_477 * DENDRO_811 + DENDRO_479 * DENDRO_810;
const double DENDRO_822 = -DENDRO_694;
const double DENDRO_823 = -DENDRO_645;
const double DENDRO_824 = DENDRO_10 * DENDRO_806;
const double DENDRO_825 =
    DENDRO_634 * DENDRO_663 + DENDRO_635 * DENDRO_811 + DENDRO_636 * DENDRO_810;
const double DENDRO_826 =
    DENDRO_634 * DENDRO_650 + DENDRO_635 * DENDRO_810 + DENDRO_636 * DENDRO_815;
const double DENDRO_827 = -DENDRO_652;
const double DENDRO_828 = -DENDRO_680;
const double DENDRO_829 =
    DENDRO_634 * DENDRO_673 + DENDRO_635 * DENDRO_663 + DENDRO_636 * DENDRO_650;
const double DENDRO_830 = -DENDRO_665;

// Dendro: printing variables
//--
psi4_real[pp] =
    (1.0 / 2.0) * DENDRO_177 * DENDRO_407 *
        (DENDRO_263 * DENDRO_658 * DENDRO_679 * chi[pp] * x +
         DENDRO_263 * DENDRO_658 * DENDRO_684 * chi[pp] * y -
         DENDRO_344 * DENDRO_657 * DENDRO_659 -
         DENDRO_344 * DENDRO_669 * DENDRO_670 +
         DENDRO_346 * DENDRO_693 * chi[pp] *
             (DENDRO_279 * DENDRO_695 + DENDRO_342 * DENDRO_697 +
              DENDRO_420 * DENDRO_654) +
         DENDRO_346 * DENDRO_701 * DENDRO_702 * chi[pp] +
         DENDRO_346 * DENDRO_705 * DENDRO_706 * chi[pp] -
         DENDRO_423 * DENDRO_687 *
             (DENDRO_279 * DENDRO_691 + DENDRO_342 * DENDRO_689 +
              DENDRO_420 * DENDRO_664) -
         DENDRO_423 * DENDRO_699 * DENDRO_700 -
         DENDRO_423 * DENDRO_703 * DENDRO_704) -
    1.0 / 8.0 * DENDRO_245 * DENDRO_280 - 1.0 / 8.0 * DENDRO_340 * DENDRO_343 +
    (1.0 / 8.0) * DENDRO_351 * DENDRO_406 +
    (1.0 / 4.0) * DENDRO_407 * DENDRO_5 *
        (DENDRO_263 * DENDRO_415 * DENDRO_417 * chi[pp] -
         DENDRO_421 * DENDRO_422 * DENDRO_423) -
    1.0 / 24.0 * DENDRO_424 * DENDRO_509 * DENDRO_510 -
    1.0 / 4.0 * DENDRO_511 * DENDRO_577 - 1.0 / 4.0 * DENDRO_578 * DENDRO_638 -
    1.0 / 24.0 * DENDRO_644 *
        (At0[pp] * DENDRO_280 + At1[pp] * DENDRO_424 * DENDRO_639 +
         At3[pp] * DENDRO_343 - At5[pp] * DENDRO_351 + DENDRO_511 * DENDRO_640 +
         DENDRO_578 * DENDRO_641);
//--
psi4_img[pp] =
    (1.0 / 12.0) *
    (6 * DENDRO_177 * DENDRO_407 *
         (-DENDRO_657 * DENDRO_701 + DENDRO_658 * DENDRO_700 * x +
          DENDRO_658 * DENDRO_704 * y - DENDRO_659 * DENDRO_702 -
          DENDRO_669 * DENDRO_705 - DENDRO_670 * DENDRO_706 +
          DENDRO_679 * DENDRO_699 + DENDRO_684 * DENDRO_703 +
          DENDRO_687 * (-DENDRO_271 * DENDRO_664 + DENDRO_277 * DENDRO_691 +
                        DENDRO_341 * DENDRO_689) -
          DENDRO_693 * (-DENDRO_271 * DENDRO_654 + DENDRO_277 * DENDRO_695 +
                        DENDRO_341 * DENDRO_697)) +
     DENDRO_233 * DENDRO_271 * DENDRO_743 *
         (-DENDRO_114 * (0.25 * DENDRO_735 + 1.0 * DENDRO_736) -
          DENDRO_114 * (-1.0 * DENDRO_223 * DENDRO_54 + DENDRO_367) +
          DENDRO_124 * (DENDRO_741 + DENDRO_742) +
          DENDRO_124 * (DENDRO_222 * DENDRO_40 + DENDRO_41 * DENDRO_46) +
          DENDRO_130 * (0.25 * DENDRO_741 + 1.0 * DENDRO_742) +
          DENDRO_130 * (-1.0 * DENDRO_192 * DENDRO_54 + DENDRO_359) -
          DENDRO_130 * (-DENDRO_222 * DENDRO_362 - DENDRO_737) -
          DENDRO_132 * (DENDRO_111 * DENDRO_219 + DENDRO_46 * DENDRO_527) -
          DENDRO_132 * (DENDRO_149 * DENDRO_207 + DENDRO_739) -
          DENDRO_132 * (DENDRO_187 * DENDRO_87 + DENDRO_740) -
          DENDRO_132 * (DENDRO_65 * DENDRO_728 + 0.25 * DENDRO_738) -
          DENDRO_178 * (DENDRO_169 * DENDRO_211 + DENDRO_173 * DENDRO_192 +
                        DENDRO_176 * DENDRO_223 + DENDRO_352) +
          DENDRO_180 * DENDRO_207 * DENDRO_59 - DENDRO_180 * DENDRO_721 +
          6.0 * DENDRO_184 * DENDRO_210 * DENDRO_50 - DENDRO_184 * DENDRO_722 +
          DENDRO_219 * DENDRO_387 - DENDRO_313 * (DENDRO_735 + DENDRO_736) -
          DENDRO_313 * (DENDRO_135 * DENDRO_222 + DENDRO_219 * DENDRO_41) +
          DENDRO_338 * DENDRO_41 * DENDRO_87 - DENDRO_382 * DENDRO_383 -
          DENDRO_384 * DENDRO_385 + DENDRO_388 * DENDRO_46 + DENDRO_405 -
          DENDRO_41 * DENDRO_723 - DENDRO_724 * gt5[pp]) -
     3 * DENDRO_277 * DENDRO_709 *
         (4 * DENDRO_0 * gt0[pp] + 4 * DENDRO_1 * gt1[pp] +
          4 * DENDRO_10 * DENDRO_137 * DENDRO_49 +
          4 * DENDRO_101 * DENDRO_33 * DENDRO_49 * DENDRO_97 - DENDRO_106 +
          4 * DENDRO_107 * DENDRO_108 * DENDRO_37 * DENDRO_49 -
          DENDRO_114 * (DENDRO_145 + DENDRO_148) -
          DENDRO_114 * (DENDRO_151 + DENDRO_153) -
          DENDRO_114 * (DENDRO_136 * DENDRO_87 + DENDRO_719) -
          DENDRO_114 * (DENDRO_46 * DENDRO_68 + DENDRO_710) -
          DENDRO_114 * (DENDRO_55 * DENDRO_64 + DENDRO_711) +
          2.0 * DENDRO_119 * DENDRO_20 * DENDRO_49 -
          DENDRO_121 * (DENDRO_715 + DENDRO_716) -
          DENDRO_121 * (DENDRO_40 * DENDRO_97 + DENDRO_720) -
          DENDRO_130 * (2 * DENDRO_108 * DENDRO_150 - DENDRO_160 * DENDRO_87) -
          DENDRO_132 * (1.0 * DENDRO_716 + DENDRO_718) -
          DENDRO_132 *
              (1.0 * DENDRO_108 * DENDRO_135 + DENDRO_160 * DENDRO_97) -
          DENDRO_141 + 4 * DENDRO_156 * DENDRO_20 * DENDRO_49 - DENDRO_159 +
          4 * DENDRO_16 * DENDRO_24 * DENDRO_25 +
          2.0 * DENDRO_16 * DENDRO_28 * DENDRO_29 +
          2.0 * DENDRO_16 * DENDRO_32 * DENDRO_33 +
          2.0 * DENDRO_16 * DENDRO_36 * DENDRO_37 -
          DENDRO_178 * (DENDRO_103 * DENDRO_173 + DENDRO_109 * DENDRO_169 +
                        DENDRO_164 + DENDRO_176 * DENDRO_71) -
          DENDRO_19 - DENDRO_23 +
          2.0 * DENDRO_25 * DENDRO_49 * (DENDRO_713 + DENDRO_714) +
          4 * DENDRO_25 * DENDRO_49 * (1.0 * DENDRO_714 + DENDRO_717) +
          2.0 * DENDRO_25 * DENDRO_49 * (DENDRO_40 * DENDRO_87 + DENDRO_712) +
          3.0 * DENDRO_29 * DENDRO_40 * DENDRO_46 * DENDRO_49 +
          4 * DENDRO_29 * DENDRO_49 * DENDRO_86 * DENDRO_88 +
          4 * DENDRO_3 * gt2[pp] +
          3.0 * DENDRO_33 * DENDRO_49 * DENDRO_52 * DENDRO_55 +
          4 * DENDRO_33 * DENDRO_49 * DENDRO_78 * DENDRO_80 +
          6.0 * DENDRO_37 * DENDRO_49 * DENDRO_60 * DENDRO_70 -
          DENDRO_40 * DENDRO_722 - DENDRO_52 * DENDRO_721 -
          DENDRO_60 * DENDRO_723 - DENDRO_724 * gt0[pp] - DENDRO_8 -
          DENDRO_96) +
     3 * DENDRO_341 * DENDRO_725 *
         (DENDRO_114 * (-DENDRO_726 - 1.0 * DENDRO_727) -
          DENDRO_114 * (DENDRO_207 * DENDRO_300 + DENDRO_302) -
          DENDRO_114 * (-1.0 * DENDRO_226 * DENDRO_45 + DENDRO_299) +
          DENDRO_121 * (DENDRO_316 - DENDRO_53 * DENDRO_55) +
          DENDRO_121 * (-DENDRO_179 * DENDRO_97 + DENDRO_315) +
          DENDRO_121 * (DENDRO_194 * DENDRO_80 - DENDRO_734) +
          DENDRO_130 * (DENDRO_289 + DENDRO_291) +
          DENDRO_130 * (DENDRO_115 * DENDRO_219 + DENDRO_294 * DENDRO_55) +
          DENDRO_130 * (DENDRO_207 * DENDRO_294 + DENDRO_731) +
          DENDRO_130 * (DENDRO_69 * DENDRO_728 + DENDRO_730) +
          DENDRO_132 * (DENDRO_308 + DENDRO_733) +
          DENDRO_132 * (0.25 * DENDRO_194 * DENDRO_80 - 1.0 * DENDRO_734) +
          DENDRO_132 * (DENDRO_214 * DENDRO_45 - DENDRO_732) -
          DENDRO_178 * (DENDRO_169 * DENDRO_214 + DENDRO_173 * DENDRO_201 +
                        DENDRO_176 * DENDRO_226 + DENDRO_281) -
          DENDRO_179 * DENDRO_722 + 6.0 * DENDRO_194 * DENDRO_200 * DENDRO_58 -
          DENDRO_194 * DENDRO_721 + DENDRO_219 * DENDRO_321 +
          DENDRO_313 * (-DENDRO_179 * DENDRO_207 + DENDRO_314) +
          DENDRO_313 * (DENDRO_182 * DENDRO_194 - DENDRO_727) +
          DENDRO_313 * (-DENDRO_219 * DENDRO_53 + DENDRO_312) +
          DENDRO_318 * DENDRO_319 - DENDRO_322 * DENDRO_323 +
          DENDRO_326 * DENDRO_55 + DENDRO_328 * DENDRO_97 + DENDRO_339 -
          DENDRO_53 * DENDRO_723 - DENDRO_724 * gt3[pp]) +
     DENDRO_5 * DENDRO_642 *
         (DENDRO_415 * DENDRO_422 + DENDRO_417 * DENDRO_421) +
     DENDRO_644 *
         (-At1[pp] * DENDRO_744 - At2[pp] * DENDRO_745 +
          At3[pp] * DENDRO_265 * DENDRO_725 * y - At4[pp] * DENDRO_746 +
          At5[pp] * DENDRO_177 * DENDRO_271 * DENDRO_743 -
          DENDRO_265 * DENDRO_408 * DENDRO_709) -
     3 * DENDRO_744 *
         (-DENDRO_114 * DENDRO_454 +
          DENDRO_114 * (-DENDRO_136 * DENDRO_200 + DENDRO_462) +
          DENDRO_114 *
              (DENDRO_298 * DENDRO_46 + DENDRO_456 - DENDRO_457 * DENDRO_55) -
          DENDRO_121 * (DENDRO_52 * DENDRO_55 + DENDRO_53 * DENDRO_70) +
          DENDRO_130 * DENDRO_442 + DENDRO_130 * (DENDRO_752 + DENDRO_753) +
          DENDRO_130 * (DENDRO_294 * DENDRO_70 + DENDRO_754) +
          DENDRO_130 * (DENDRO_749 + DENDRO_750 + DENDRO_751) +
          DENDRO_132 * (-DENDRO_200 * DENDRO_68 + DENDRO_471) -
          DENDRO_132 * (DENDRO_108 * DENDRO_198 + DENDRO_467 * DENDRO_97 +
                        DENDRO_502 * DENDRO_97) +
          DENDRO_132 * (-DENDRO_108 * DENDRO_199 +
                        0.5 * DENDRO_214 * DENDRO_40 - DENDRO_457 * DENDRO_97) +
          DENDRO_132 *
              (DENDRO_144 * DENDRO_70 - DENDRO_157 * DENDRO_55 + DENDRO_465) -
          DENDRO_178 * (DENDRO_154 * DENDRO_173 + DENDRO_176 * DENDRO_297 +
                        DENDRO_474 + DENDRO_476 * DENDRO_98) +
          DENDRO_325 * (DENDRO_108 * DENDRO_294 + DENDRO_111 * DENDRO_97 +
                        0.5 * DENDRO_720) +
          DENDRO_325 *
              (DENDRO_68 * DENDRO_70 + DENDRO_69 * DENDRO_70 + DENDRO_718) -
          DENDRO_43 * DENDRO_722 + DENDRO_433 * DENDRO_89 +
          DENDRO_435 * DENDRO_81 + DENDRO_481 + DENDRO_483 + DENDRO_484 +
          DENDRO_486 + DENDRO_487 + DENDRO_489 + DENDRO_490 + DENDRO_491 +
          DENDRO_492 + DENDRO_493 + DENDRO_494 + DENDRO_495 + DENDRO_496 +
          DENDRO_50 * (DENDRO_738 + DENDRO_747 + DENDRO_748) + DENDRO_508 -
          DENDRO_66 * DENDRO_723 - DENDRO_721 * DENDRO_75 -
          DENDRO_724 * gt1[pp] -
          DENDRO_756 * (DENDRO_46 * DENDRO_53 + DENDRO_757) -
          DENDRO_756 *
              (DENDRO_179 * DENDRO_87 + DENDRO_180 * DENDRO_97 + DENDRO_755) -
          DENDRO_81 * (DENDRO_438 - DENDRO_55 * DENDRO_74 + DENDRO_733) -
          DENDRO_81 * (0.5 * DENDRO_135 * DENDRO_214 - DENDRO_333 * DENDRO_97 -
                       DENDRO_732)) -
     3 * DENDRO_745 *
         (2.0 * DENDRO_0 * gt2[pp] + 2.0 * DENDRO_1 * gt4[pp] +
          DENDRO_10 * DENDRO_49 * DENDRO_514 +
          4 * DENDRO_10 * DENDRO_49 *
              (0.5 * DENDRO_150 * DENDRO_207 - DENDRO_760) +
          4 * DENDRO_10 * DENDRO_49 *
              (-DENDRO_222 * DENDRO_68 - DENDRO_46 * DENDRO_467 - DENDRO_764) -
          DENDRO_114 * DENDRO_542 -
          DENDRO_114 * (DENDRO_136 * DENDRO_210 + DENDRO_760) -
          DENDRO_130 * (0.5 * DENDRO_150 * DENDRO_87 - DENDRO_762) -
          DENDRO_130 * (DENDRO_150 * DENDRO_70 - DENDRO_160 * DENDRO_46 -
                        DENDRO_222 * DENDRO_61) -
          DENDRO_132 * (DENDRO_749 + DENDRO_758) -
          DENDRO_132 * (DENDRO_753 + DENDRO_765) -
          DENDRO_132 * (DENDRO_502 * DENDRO_87 + DENDRO_758) -
          DENDRO_132 * (DENDRO_527 * DENDRO_70 + DENDRO_754) +
          4 * DENDRO_16 * DENDRO_25 * DENDRO_570 +
          2.0 * DENDRO_16 * DENDRO_29 * DENDRO_571 +
          2.0 * DENDRO_16 * DENDRO_33 * DENDRO_572 +
          2.0 * DENDRO_16 * DENDRO_37 * DENDRO_573 -
          DENDRO_178 * (DENDRO_176 * DENDRO_360 + DENDRO_476 * DENDRO_88 +
                        DENDRO_478 * DENDRO_91 + DENDRO_561) +
          4 * DENDRO_20 * DENDRO_49 * DENDRO_551 +
          4 * DENDRO_20 * DENDRO_49 * DENDRO_556 +
          4 * DENDRO_25 * DENDRO_49 * (DENDRO_210 * DENDRO_64 + DENDRO_762) +
          2.0 * DENDRO_25 * DENDRO_49 *
              (DENDRO_40 * DENDRO_46 + DENDRO_41 * DENDRO_70) +
          4 * DENDRO_29 * DENDRO_49 * DENDRO_518 + 2.0 * DENDRO_3 * gt5[pp] +
          DENDRO_33 * DENDRO_49 * (DENDRO_757 + DENDRO_763) +
          4 * DENDRO_33 * DENDRO_49 *
              (DENDRO_187 * DENDRO_97 + 0.25 * DENDRO_755) +
          4 * DENDRO_37 * DENDRO_49 *
              (DENDRO_111 * DENDRO_87 + 1.0 * DENDRO_712) +
          4 * DENDRO_37 * DENDRO_49 *
              (DENDRO_64 * DENDRO_70 + DENDRO_65 * DENDRO_70 + DENDRO_717) +
          2.0 * DENDRO_391 * gt0[pp] + 2.0 * DENDRO_392 * gt1[pp] +
          2.0 * DENDRO_393 * gt2[pp] - DENDRO_44 * DENDRO_721 - DENDRO_520 -
          DENDRO_524 - DENDRO_529 - DENDRO_533 - DENDRO_536 - DENDRO_569 -
          DENDRO_574 - DENDRO_575 - DENDRO_62 * DENDRO_723 -
          DENDRO_722 * DENDRO_83 - DENDRO_724 * gt2[pp] -
          DENDRO_756 * (DENDRO_738 + DENDRO_761) -
          DENDRO_89 *
              (-DENDRO_222 * DENDRO_64 - DENDRO_46 * DENDRO_82 - DENDRO_737)) -
     3 * DENDRO_746 *
         (DENDRO_114 * (DENDRO_191 * DENDRO_200 + DENDRO_596) -
          DENDRO_114 *
              (DENDRO_143 * DENDRO_222 + DENDRO_219 * DENDRO_467 + DENDRO_593) +
          DENDRO_114 * (0.5 * DENDRO_184 * DENDRO_214 -
                        DENDRO_198 * DENDRO_210 - DENDRO_767) +
          DENDRO_114 * (DENDRO_207 * DENDRO_358 + DENDRO_599 - DENDRO_767) +
          DENDRO_114 * (-DENDRO_219 * DENDRO_457 - DENDRO_219 * DENDRO_502 +
                        0.5 * DENDRO_226 * DENDRO_41) +
          DENDRO_123 * (DENDRO_748 + DENDRO_761) -
          DENDRO_130 * (0.5 * DENDRO_191 * DENDRO_87 - DENDRO_766) +
          DENDRO_130 * (DENDRO_210 * DENDRO_294 + DENDRO_766) -
          DENDRO_130 *
              (-DENDRO_160 * DENDRO_219 - DENDRO_222 * DENDRO_69 - DENDRO_764) +
          DENDRO_132 * (-DENDRO_200 * DENDRO_294 + DENDRO_612) +
          DENDRO_132 * (-DENDRO_200 * DENDRO_527 + DENDRO_610) +
          DENDRO_132 * (-DENDRO_207 * DENDRO_457 - DENDRO_768) +
          DENDRO_132 * (-DENDRO_333 * DENDRO_87 - DENDRO_768) +
          DENDRO_132 *
              (-DENDRO_157 * DENDRO_219 + DENDRO_226 * DENDRO_65 - DENDRO_730) +
          DENDRO_132 * (0.5 * DENDRO_226 * DENDRO_40 - DENDRO_46 * DENDRO_74 -
                        0.25 * DENDRO_763) -
          DENDRO_178 * (DENDRO_173 * DENDRO_400 + DENDRO_208 * DENDRO_476 +
                        DENDRO_220 * DENDRO_480 + DENDRO_633) -
          DENDRO_188 * DENDRO_722 - DENDRO_196 * DENDRO_721 +
          DENDRO_313 * (DENDRO_179 * DENDRO_182 - DENDRO_180 * DENDRO_200) +
          DENDRO_325 * (DENDRO_149 * DENDRO_97 + DENDRO_751) +
          DENDRO_325 * (DENDRO_710 + DENDRO_752 + DENDRO_765) -
          DENDRO_42 * DENDRO_723 + DENDRO_580 * DENDRO_89 + DENDRO_632 -
          DENDRO_724 * gt4[pp] -
          DENDRO_81 * (-DENDRO_207 * DENDRO_333 + DENDRO_583) -
          DENDRO_81 * (-DENDRO_219 * DENDRO_74 + DENDRO_587) -
          DENDRO_81 * (-DENDRO_198 * DENDRO_200 - DENDRO_199 * DENDRO_200 -
                       DENDRO_726) -
          DENDRO_89 * (0.5 * DENDRO_150 * DENDRO_219 - DENDRO_219 * DENDRO_82 -
                       DENDRO_222 * DENDRO_294))) /
    (sqrt(DENDRO_177 * DENDRO_262) * sqrt(DENDRO_177 * DENDRO_345));
//--
ham[pp] =
    -At0[pp] * DENDRO_49 *
        (At0[pp] * (DENDRO_37 * DENDRO_37) + At3[pp] * DENDRO_769 +
         At5[pp] * DENDRO_770 - DENDRO_239 * DENDRO_771 +
         DENDRO_239 * DENDRO_774 - DENDRO_772 * DENDRO_773) -
    2 * At1[pp] * DENDRO_49 *
        (At1[pp] * DENDRO_33 * DENDRO_37 + At1[pp] * DENDRO_769 -
         At2[pp] * DENDRO_782 + At4[pp] * DENDRO_10 * DENDRO_20 +
         At4[pp] * DENDRO_25 * DENDRO_33 - DENDRO_10 * DENDRO_774 -
         DENDRO_10 * DENDRO_783 - DENDRO_20 * DENDRO_781 -
         DENDRO_20 * DENDRO_784) -
    At3[pp] * DENDRO_49 *
        (At0[pp] * DENDRO_769 - At1[pp] * DENDRO_33 * DENDRO_772 +
         At3[pp] * (DENDRO_33 * DENDRO_33) - At4[pp] * DENDRO_33 * DENDRO_777 +
         At5[pp] * DENDRO_775 + DENDRO_772 * DENDRO_776) -
    At5[pp] * DENDRO_49 *
        (At0[pp] * DENDRO_770 + At3[pp] * DENDRO_775 +
         At5[pp] * (DENDRO_29 * DENDRO_29) - DENDRO_239 * DENDRO_778 +
         DENDRO_239 * DENDRO_779 - DENDRO_777 * DENDRO_780) +
    (1.0 / 4.0) * DENDRO_10 * DENDRO_16 * DENDRO_638 * chi[pp] +
    (1.0 / 4.0) * DENDRO_10 * DENDRO_16 * chi[pp] *
        (DENDRO_114 * (DENDRO_594 + DENDRO_597) +
         DENDRO_114 * (DENDRO_596 + DENDRO_800) +
         DENDRO_114 * (-DENDRO_144 * DENDRO_223 + DENDRO_598) +
         DENDRO_114 * (DENDRO_187 * DENDRO_201 + DENDRO_595 + DENDRO_800) +
         DENDRO_114 * (DENDRO_199 * DENDRO_211 + DENDRO_599 + DENDRO_601) +
         DENDRO_120 * (DENDRO_290 + DENDRO_292 + DENDRO_425) +
         DENDRO_120 * (DENDRO_295 + DENDRO_428 + DENDRO_559) -
         DENDRO_130 * (DENDRO_373 + DENDRO_589) -
         DENDRO_130 * (DENDRO_399 + DENDRO_623) -
         DENDRO_130 * (DENDRO_399 + DENDRO_792) -
         DENDRO_130 * (DENDRO_590 + DENDRO_793) -
         DENDRO_130 * (DENDRO_592 + DENDRO_794) -
         DENDRO_130 * (DENDRO_795 + DENDRO_799) +
         DENDRO_132 * (DENDRO_335 + DENDRO_612) +
         DENDRO_132 * (DENDRO_455 + DENDRO_602) +
         DENDRO_132 * (DENDRO_335 + DENDRO_461 + DENDRO_609) +
         DENDRO_132 * (DENDRO_300 * DENDRO_88 + DENDRO_603 + DENDRO_604) +
         DENDRO_313 * (DENDRO_179 * DENDRO_211 + DENDRO_389) -
         DENDRO_325 * (DENDRO_116 + DENDRO_588) -
         DENDRO_325 * (DENDRO_143 * DENDRO_91 + DENDRO_555) + DENDRO_614 +
         DENDRO_615 + DENDRO_616 + DENDRO_617 + DENDRO_618 + DENDRO_619 +
         DENDRO_620 + DENDRO_624 + DENDRO_625 + DENDRO_626 + DENDRO_627 -
         DENDRO_629 * (DENDRO_445 + DENDRO_789) + DENDRO_630 + DENDRO_631 +
         DENDRO_637 - DENDRO_81 * (0.5 * DENDRO_303 + DENDRO_584) -
         DENDRO_81 * (-DENDRO_299 + DENDRO_585 + DENDRO_586) -
         DENDRO_81 * (DENDRO_187 * DENDRO_214 + DENDRO_301 + DENDRO_582) -
         DENDRO_89 * (0.5 * DENDRO_379 + DENDRO_581) -
         DENDRO_89 * (1.0 * DENDRO_402 + DENDRO_621) -
         DENDRO_89 * (DENDRO_187 * DENDRO_211 + DENDRO_369 - DENDRO_579)) +
    (1.0 / 4.0) * DENDRO_16 * DENDRO_20 * DENDRO_509 * chi[pp] +
    (1.0 / 4.0) * DENDRO_16 * DENDRO_20 * chi[pp] *
        (DENDRO_114 * (-DENDRO_291 - DENDRO_453) +
         DENDRO_114 * (DENDRO_296 + DENDRO_458) +
         DENDRO_114 * (DENDRO_296 + DENDRO_607) +
         DENDRO_114 * (DENDRO_459 + DENDRO_610) +
         DENDRO_114 * (DENDRO_461 + DENDRO_505 + DENDRO_611) +
         DENDRO_114 * (DENDRO_522 + DENDRO_604 + DENDRO_608) +
         DENDRO_121 * (DENDRO_201 * DENDRO_52 + DENDRO_337) -
         DENDRO_130 * (DENDRO_447 + DENDRO_801) -
         DENDRO_130 * (DENDRO_449 + DENDRO_801) -
         DENDRO_130 * (DENDRO_503 + DENDRO_796) -
         DENDRO_130 * (0.5 * DENDRO_208 * DENDRO_65 - DENDRO_441) +
         DENDRO_132 * (DENDRO_464 + DENDRO_468) +
         DENDRO_132 * (DENDRO_466 + DENDRO_802) +
         DENDRO_132 * (DENDRO_214 * DENDRO_65 + DENDRO_469) +
         DENDRO_132 * (DENDRO_143 * DENDRO_71 + DENDRO_465 + DENDRO_802) +
         DENDRO_132 * (DENDRO_201 * DENDRO_69 + DENDRO_470 + DENDRO_507) -
         DENDRO_325 * (0.5 * DENDRO_127 + DENDRO_439) -
         DENDRO_325 * (DENDRO_162 + DENDRO_440) -
         DENDRO_325 * (DENDRO_103 * DENDRO_143 + DENDRO_155 + DENDRO_500) -
         DENDRO_473 * (DENDRO_396 + DENDRO_398 + DENDRO_498) + DENDRO_497 -
         DENDRO_613 * (DENDRO_133 + DENDRO_146 + DENDRO_790) -
         DENDRO_613 * (DENDRO_445 + DENDRO_787 + DENDRO_788) -
         DENDRO_81 * (0.5 * DENDRO_315 + DENDRO_436) -
         DENDRO_81 * (1.0 * DENDRO_316 + DENDRO_437) -
         DENDRO_81 * (DENDRO_143 * DENDRO_201 + DENDRO_310 - DENDRO_434) -
         DENDRO_89 * (0.5 * DENDRO_376 + DENDRO_799) -
         DENDRO_89 * (0.25 * DENDRO_208 * DENDRO_41 - DENDRO_432)) -
    DENDRO_245 * DENDRO_38 * DENDRO_786 - DENDRO_30 * DENDRO_406 * DENDRO_786 -
    DENDRO_34 * DENDRO_340 * DENDRO_786 -
    DENDRO_49 * DENDRO_640 *
        (-At1[pp] * DENDRO_782 + At2[pp] * DENDRO_770 +
         At3[pp] * DENDRO_10 * DENDRO_20 - At4[pp] * DENDRO_10 * DENDRO_25 -
         DENDRO_10 * DENDRO_773 - DENDRO_20 * DENDRO_780 +
         DENDRO_25 * DENDRO_781 + DENDRO_29 * DENDRO_774 +
         DENDRO_29 * DENDRO_783) -
    DENDRO_49 * DENDRO_641 *
        (-At0[pp] * DENDRO_782 + At1[pp] * DENDRO_10 * DENDRO_20 +
         At1[pp] * DENDRO_25 * DENDRO_33 + At4[pp] * DENDRO_29 * DENDRO_33 +
         At4[pp] * DENDRO_775 - DENDRO_10 * DENDRO_784 -
         DENDRO_20 * DENDRO_779 - DENDRO_25 * DENDRO_776 -
         DENDRO_29 * DENDRO_785) -
    DENDRO_577 * DENDRO_798 -
    DENDRO_798 *
        (DENDRO_114 * (DENDRO_377 + DENDRO_545) +
         DENDRO_114 * (DENDRO_377 + DENDRO_795) +
         DENDRO_114 * (DENDRO_499 + DENDRO_792) +
         DENDRO_114 * (DENDRO_546 + DENDRO_793) +
         DENDRO_114 * (DENDRO_549 + DENDRO_794) +
         DENDRO_114 * (0.5 * DENDRO_395 + DENDRO_397 - DENDRO_541) +
         DENDRO_120 * (DENDRO_628 + DENDRO_790) +
         DENDRO_120 * (DENDRO_787 + DENDRO_789) -
         DENDRO_124 * (DENDRO_211 * DENDRO_40 + DENDRO_390) -
         DENDRO_130 * (DENDRO_531 + DENDRO_791) -
         DENDRO_130 * (DENDRO_532 + DENDRO_534) -
         DENDRO_130 * (DENDRO_192 * DENDRO_69 + DENDRO_535) -
         DENDRO_130 * (DENDRO_149 * DENDRO_71 + DENDRO_530 + DENDRO_791) -
         DENDRO_130 * (DENDRO_211 * DENDRO_65 + DENDRO_537 + DENDRO_539) +
         DENDRO_132 * (DENDRO_550 + DENDRO_796) +
         DENDRO_132 * (DENDRO_557 + DENDRO_797) +
         DENDRO_132 * (DENDRO_558 + DENDRO_797) +
         DENDRO_132 * (DENDRO_160 * DENDRO_208 + DENDRO_444 + DENDRO_552) -
         DENDRO_325 * (0.5 * DENDRO_125 + DENDRO_526) -
         DENDRO_325 * (DENDRO_158 + DENDRO_528) -
         DENDRO_325 * (DENDRO_109 * DENDRO_149 + DENDRO_161 + DENDRO_525) -
         DENDRO_560 * (DENDRO_290 + DENDRO_427) + DENDRO_576 -
         DENDRO_81 * (DENDRO_220 * DENDRO_68 + DENDRO_606) -
         DENDRO_81 * (DENDRO_335 + DENDRO_505 + DENDRO_523) -
         DENDRO_89 * (1.0 * DENDRO_381 + DENDRO_521) -
         DENDRO_89 * (0.5 * DENDRO_404 + DENDRO_519) -
         DENDRO_89 * (DENDRO_149 * DENDRO_211 + DENDRO_364 - DENDRO_517)) +
    (2.0 / 3.0) * pow(K[pp], 2);
//--
mom0[pp] =
    (3.0 / 2.0) * At0[pp] * DENDRO_16 * DENDRO_168 * DENDRO_177 * DENDRO_25 +
    (3.0 / 2.0) * At0[pp] * DENDRO_16 * DENDRO_177 * DENDRO_37 * DENDRO_6 -
    At0[pp] * DENDRO_805 - At0[pp] * Gt0[pp] +
    (3.0 / 2.0) * At1[pp] * DENDRO_16 * DENDRO_172 * DENDRO_177 * DENDRO_33 -
    At1[pp] * DENDRO_807 - At1[pp] * Gt1[pp] +
    (3.0 / 2.0) * At2[pp] * DENDRO_16 * DENDRO_168 * DENDRO_177 * DENDRO_29 +
    (3.0 / 2.0) * At2[pp] * DENDRO_16 * DENDRO_177 * DENDRO_25 * DENDRO_6 -
    At2[pp] * Gt2[pp] + DENDRO_16 * DENDRO_25 * (-DENDRO_677 + DENDRO_817) +
    DENDRO_16 * DENDRO_25 *
        (DENDRO_167 * DENDRO_672 + DENDRO_171 * DENDRO_662 +
         DENDRO_175 * DENDRO_649 + DENDRO_823) +
    DENDRO_16 * DENDRO_29 * (DENDRO_818 + DENDRO_819) +
    DENDRO_16 * DENDRO_33 * (DENDRO_820 + DENDRO_821) +
    DENDRO_16 * DENDRO_37 *
        (DENDRO_167 * DENDRO_649 + DENDRO_171 * DENDRO_646 +
         DENDRO_175 * DENDRO_648 - grad_0_At0[pp]) -
    DENDRO_17 * (DENDRO_809 + DENDRO_812) -
    DENDRO_17 * (DENDRO_813 + DENDRO_814) -
    DENDRO_21 * (-DENDRO_690 + DENDRO_816) -
    DENDRO_21 * (DENDRO_167 * DENDRO_662 + DENDRO_171 * DENDRO_661 +
                 DENDRO_175 * DENDRO_646 + DENDRO_822) -
    DENDRO_776 * DENDRO_808 - DENDRO_778 * DENDRO_806 -
    2.0 / 3.0 * grad_0_K[pp];
//--
mom1[pp] =
    (3.0 / 2.0) * At1[pp] * DENDRO_16 * DENDRO_168 * DENDRO_177 * DENDRO_25 +
    (3.0 / 2.0) * At1[pp] * DENDRO_16 * DENDRO_177 * DENDRO_37 * DENDRO_6 -
    At1[pp] * DENDRO_805 - At1[pp] * Gt0[pp] +
    (3.0 / 2.0) * At3[pp] * DENDRO_16 * DENDRO_172 * DENDRO_177 * DENDRO_33 -
    At3[pp] * DENDRO_807 - At3[pp] * DENDRO_824 - At3[pp] * Gt1[pp] -
    At4[pp] * DENDRO_10 * DENDRO_808 +
    (3.0 / 2.0) * At4[pp] * DENDRO_16 * DENDRO_168 * DENDRO_177 * DENDRO_29 +
    (3.0 / 2.0) * At4[pp] * DENDRO_16 * DENDRO_177 * DENDRO_25 * DENDRO_6 -
    At4[pp] * Gt2[pp] + DENDRO_16 * DENDRO_25 * (DENDRO_809 + DENDRO_826) +
    DENDRO_16 * DENDRO_25 * (DENDRO_814 + DENDRO_827) +
    DENDRO_16 * DENDRO_29 * (DENDRO_828 + DENDRO_829) +
    DENDRO_16 * DENDRO_33 *
        (DENDRO_283 * DENDRO_662 + DENDRO_284 * DENDRO_661 +
         DENDRO_285 * DENDRO_646 - grad_1_At3[pp]) +
    DENDRO_16 * DENDRO_37 * (DENDRO_816 + DENDRO_822) -
    DENDRO_17 * (-DENDRO_682 + DENDRO_825) -
    DENDRO_17 * (DENDRO_283 * DENDRO_672 + DENDRO_284 * DENDRO_662 +
                 DENDRO_285 * DENDRO_649 + DENDRO_830) -
    DENDRO_21 * (-DENDRO_696 + DENDRO_821) -
    DENDRO_21 * (DENDRO_283 * DENDRO_649 + DENDRO_284 * DENDRO_646 +
                 DENDRO_285 * DENDRO_648 + DENDRO_820) -
    2.0 / 3.0 * grad_1_K[pp];
//--
mom2[pp] =
    (3.0 / 2.0) * At2[pp] * DENDRO_16 * DENDRO_168 * DENDRO_177 * DENDRO_25 +
    (3.0 / 2.0) * At2[pp] * DENDRO_16 * DENDRO_177 * DENDRO_37 * DENDRO_6 -
    At2[pp] * DENDRO_805 - At2[pp] * Gt0[pp] +
    (3.0 / 2.0) * At4[pp] * DENDRO_16 * DENDRO_172 * DENDRO_177 * DENDRO_33 -
    At4[pp] * DENDRO_824 - At4[pp] * Gt1[pp] +
    (3.0 / 2.0) * At5[pp] * DENDRO_16 * DENDRO_168 * DENDRO_177 * DENDRO_29 +
    (3.0 / 2.0) * At5[pp] * DENDRO_16 * DENDRO_177 * DENDRO_25 * DENDRO_6 -
    At5[pp] * Gt2[pp] + DENDRO_16 * DENDRO_25 * (-DENDRO_655 + DENDRO_819) +
    DENDRO_16 * DENDRO_25 *
        (DENDRO_355 * DENDRO_649 + DENDRO_356 * DENDRO_646 +
         DENDRO_357 * DENDRO_648 + DENDRO_818) +
    DENDRO_16 * DENDRO_29 *
        (DENDRO_355 * DENDRO_672 + DENDRO_356 * DENDRO_662 +
         DENDRO_357 * DENDRO_649 - grad_2_At5[pp]) +
    DENDRO_16 * DENDRO_33 * (DENDRO_825 + DENDRO_830) +
    DENDRO_16 * DENDRO_37 * (DENDRO_817 + DENDRO_823) -
    DENDRO_17 * (-DENDRO_667 + DENDRO_829) -
    DENDRO_17 * (DENDRO_355 * DENDRO_662 + DENDRO_356 * DENDRO_661 +
                 DENDRO_357 * DENDRO_646 + DENDRO_828) -
    DENDRO_176 * DENDRO_771 * DENDRO_803 -
    DENDRO_21 * (DENDRO_812 + DENDRO_827) -
    DENDRO_21 * (DENDRO_813 + DENDRO_826) - DENDRO_785 * DENDRO_808 -
    2.0 / 3.0 * grad_2_K[pp];
// Dendro: reduced ops: 3968
// Dendro: }}}
