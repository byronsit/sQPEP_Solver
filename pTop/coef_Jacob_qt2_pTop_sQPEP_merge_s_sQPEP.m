function obj = coef_Jacob_qt2_pTop_sQPEP_merge_s_sQPEP(in1,in2,in3)
coef_J2 = in3(:,2);
coef_J5 = in3(:,5);
coef_J6 = in3(:,6);
coef_J7 = in3(:,7);
coef_J15 = in3(:,15);
coef_J16 = in3(:,16);
coef_J17 = in3(:,17);
coef_J18 = in3(:,18);
coef_J19 = in3(:,19);
coef_J20 = in3(:,20);
coef_J21 = in3(:,21);
coef_J22 = in3(:,22);
coef_J23 = in3(:,23);
coef_J24 = in3(:,24);
coef_J37 = in3(:,37);
coef_J38 = in3(:,38);
coef_J39 = in3(:,39);
coef_J40 = in3(:,40);
coef_J41 = in3(:,41);
coef_J42 = in3(:,42);
coef_J43 = in3(:,43);
coef_J44 = in3(:,44);
coef_J45 = in3(:,45);
coef_J46 = in3(:,46);
coef_J47 = in3(:,47);
coef_J48 = in3(:,48);
coef_J49 = in3(:,49);
coef_J50 = in3(:,50);
coef_J51 = in3(:,51);
coef_J52 = in3(:,52);
coef_J53 = in3(:,53);
coef_J54 = in3(:,54);
coef_J55 = in3(:,55);
coef_J56 = in3(:,56);
coef_J57 = in3(:,57);
coef_J58 = in3(:,58);
coefs_tq1_1 = in1(1);
coefs_tq1_2 = in1(4);
coefs_tq1_3 = in1(7);
coefs_tq1_4 = in1(10);
coefs_tq1_5 = in1(13);
coefs_tq1_6 = in1(16);
coefs_tq1_7 = in1(19);
coefs_tq1_8 = in1(22);
coefs_tq1_9 = in1(25);
coefs_tq2_1 = in1(2);
coefs_tq2_2 = in1(5);
coefs_tq2_3 = in1(8);
coefs_tq2_4 = in1(11);
coefs_tq2_5 = in1(14);
coefs_tq2_6 = in1(17);
coefs_tq2_7 = in1(20);
coefs_tq2_8 = in1(23);
coefs_tq2_9 = in1(26);
coefs_tq3_1 = in1(3);
coefs_tq3_2 = in1(6);
coefs_tq3_3 = in1(9);
coefs_tq3_4 = in1(12);
coefs_tq3_5 = in1(15);
coefs_tq3_6 = in1(18);
coefs_tq3_7 = in1(21);
coefs_tq3_8 = in1(24);
coefs_tq3_9 = in1(27);
coefs_tq1_10 = in1(28);
coefs_tq1_11 = in1(31);
coefs_tq2_10 = in1(29);
coefs_tq2_11 = in1(32);
coefs_tq3_10 = in1(30);
coefs_tq3_11 = in1(33);
pinvG1_1 = in2(1);
pinvG1_2 = in2(4);
pinvG1_3 = in2(7);
pinvG2_1 = in2(2);
pinvG2_2 = in2(5);
pinvG2_3 = in2(8);
pinvG3_1 = in2(3);
pinvG3_2 = in2(6);
pinvG3_3 = in2(9);
t2 = coefs_tq1_1.*pinvG1_1;
t3 = coefs_tq1_2.*pinvG1_1;
t4 = coefs_tq1_3.*pinvG1_1;
t5 = coefs_tq1_4.*pinvG1_1;
t6 = coefs_tq1_5.*pinvG1_1;
t7 = coefs_tq1_6.*pinvG1_1;
t8 = coefs_tq1_7.*pinvG1_1;
t9 = coefs_tq1_8.*pinvG1_1;
t10 = coefs_tq1_9.*pinvG1_1;
t11 = coefs_tq1_1.*pinvG2_1;
t12 = coefs_tq1_2.*pinvG2_1;
t13 = coefs_tq2_1.*pinvG1_2;
t14 = coefs_tq1_3.*pinvG2_1;
t15 = coefs_tq2_2.*pinvG1_2;
t16 = coefs_tq1_4.*pinvG2_1;
t17 = coefs_tq2_3.*pinvG1_2;
t18 = coefs_tq1_5.*pinvG2_1;
t19 = coefs_tq2_4.*pinvG1_2;
t20 = coefs_tq1_6.*pinvG2_1;
t21 = coefs_tq2_5.*pinvG1_2;
t22 = coefs_tq1_7.*pinvG2_1;
t23 = coefs_tq2_6.*pinvG1_2;
t24 = coefs_tq1_8.*pinvG2_1;
t25 = coefs_tq2_7.*pinvG1_2;
t26 = coefs_tq1_9.*pinvG2_1;
t27 = coefs_tq2_8.*pinvG1_2;
t28 = coefs_tq2_9.*pinvG1_2;
t29 = coefs_tq1_1.*pinvG3_1;
t30 = coefs_tq1_2.*pinvG3_1;
t31 = coefs_tq2_1.*pinvG2_2;
t32 = coefs_tq1_3.*pinvG3_1;
t33 = coefs_tq2_2.*pinvG2_2;
t34 = coefs_tq3_1.*pinvG1_3;
t35 = coefs_tq1_4.*pinvG3_1;
t36 = coefs_tq2_3.*pinvG2_2;
t37 = coefs_tq3_2.*pinvG1_3;
t38 = coefs_tq1_5.*pinvG3_1;
t39 = coefs_tq2_4.*pinvG2_2;
t40 = coefs_tq3_3.*pinvG1_3;
t41 = coefs_tq1_6.*pinvG3_1;
t42 = coefs_tq2_5.*pinvG2_2;
t43 = coefs_tq3_4.*pinvG1_3;
t44 = coefs_tq1_7.*pinvG3_1;
t45 = coefs_tq2_6.*pinvG2_2;
t46 = coefs_tq3_5.*pinvG1_3;
t47 = coefs_tq1_8.*pinvG3_1;
t48 = coefs_tq2_7.*pinvG2_2;
t49 = coefs_tq3_6.*pinvG1_3;
t50 = coefs_tq1_9.*pinvG3_1;
t51 = coefs_tq2_8.*pinvG2_2;
t52 = coefs_tq3_7.*pinvG1_3;
t53 = coefs_tq2_9.*pinvG2_2;
t54 = coefs_tq3_8.*pinvG1_3;
t55 = coefs_tq3_9.*pinvG1_3;
t56 = coefs_tq2_1.*pinvG3_2;
t57 = coefs_tq2_2.*pinvG3_2;
t58 = coefs_tq3_1.*pinvG2_3;
t59 = coefs_tq2_3.*pinvG3_2;
t60 = coefs_tq3_2.*pinvG2_3;
t61 = coefs_tq2_4.*pinvG3_2;
t62 = coefs_tq3_3.*pinvG2_3;
t63 = coefs_tq2_5.*pinvG3_2;
t64 = coefs_tq3_4.*pinvG2_3;
t65 = coefs_tq2_6.*pinvG3_2;
t66 = coefs_tq3_5.*pinvG2_3;
t67 = coefs_tq2_7.*pinvG3_2;
t68 = coefs_tq3_6.*pinvG2_3;
t69 = coefs_tq2_8.*pinvG3_2;
t70 = coefs_tq3_7.*pinvG2_3;
t71 = coefs_tq2_9.*pinvG3_2;
t72 = coefs_tq3_8.*pinvG2_3;
t73 = coefs_tq3_9.*pinvG2_3;
t74 = coefs_tq3_1.*pinvG3_3;
t75 = coefs_tq3_2.*pinvG3_3;
t76 = coefs_tq3_3.*pinvG3_3;
t77 = coefs_tq3_4.*pinvG3_3;
t78 = coefs_tq3_5.*pinvG3_3;
t79 = coefs_tq3_6.*pinvG3_3;
t80 = coefs_tq3_7.*pinvG3_3;
t81 = coefs_tq3_8.*pinvG3_3;
t82 = coefs_tq3_9.*pinvG3_3;
t83 = coefs_tq1_10.*pinvG1_1;
t84 = coefs_tq1_11.*pinvG1_1;
t85 = coefs_tq1_10.*pinvG2_1;
t86 = coefs_tq1_11.*pinvG2_1;
t87 = coefs_tq1_10.*pinvG3_1;
t88 = coefs_tq1_11.*pinvG3_1;
t89 = coefs_tq2_10.*pinvG1_2;
t90 = coefs_tq2_11.*pinvG1_2;
t91 = coefs_tq2_10.*pinvG2_2;
t92 = coefs_tq2_11.*pinvG2_2;
t93 = coefs_tq2_10.*pinvG3_2;
t94 = coefs_tq2_11.*pinvG3_2;
t95 = coefs_tq3_10.*pinvG1_3;
t96 = coefs_tq3_11.*pinvG1_3;
t97 = coefs_tq3_10.*pinvG2_3;
t98 = coefs_tq3_11.*pinvG2_3;
t99 = coefs_tq3_10.*pinvG3_3;
t100 = coefs_tq3_11.*pinvG3_3;
t101 = t2+t13+t34;
t102 = t3+t15+t37;
t103 = t4+t17+t40;
t104 = t5+t19+t43;
t105 = t6+t21+t46;
t106 = t7+t23+t49;
t107 = t8+t25+t52;
t108 = t9+t27+t54;
t109 = t10+t28+t55;
t110 = t11+t31+t58;
t111 = t12+t33+t60;
t112 = t14+t36+t62;
t113 = t16+t39+t64;
t114 = t18+t42+t66;
t115 = t20+t45+t68;
t116 = t22+t48+t70;
t117 = t24+t51+t72;
t118 = t26+t53+t73;
t119 = t29+t56+t74;
t120 = t30+t57+t75;
t121 = t32+t59+t76;
t122 = t35+t61+t77;
t123 = t38+t63+t78;
t124 = t41+t65+t79;
t125 = t44+t67+t80;
t126 = t47+t69+t81;
t127 = t50+t71+t82;
t128 = t83+t89+t95;
t129 = t84+t90+t96;
t130 = t85+t91+t97;
t131 = t86+t92+t98;
t132 = t87+t93+t99;
t133 = t88+t94+t100;
mt1 = [coef_J2+coef_J21.*t101+coef_J22.*t110+coef_J23.*t119,coef_J5.*2.0+coef_J21.*t102+coef_J22.*t111+coef_J23.*t120+coef_J43.*t101.*2.0+coef_J44.*t110.*2.0+coef_J45.*t119.*2.0,coef_J6+coef_J21.*t103+coef_J22.*t112+coef_J23.*t121+coef_J50.*t101+coef_J51.*t110+coef_J52.*t119,coef_J7+coef_J21.*t104+coef_J22.*t113+coef_J23.*t122+coef_J55.*t101+coef_J56.*t110+coef_J57.*t119,coef_J15.*3.0+coef_J21.*t105+coef_J22.*t114+coef_J43.*t102.*2.0+coef_J23.*t123+coef_J44.*t111.*2.0+coef_J45.*t120.*2.0,coef_J16.*2.0+coef_J21.*t106+coef_J22.*t115+coef_J43.*t103.*2.0+coef_J23.*t124+coef_J50.*t102+coef_J44.*t112.*2.0+coef_J51.*t111+coef_J45.*t121.*2.0+coef_J52.*t120,coef_J17.*2.0+coef_J21.*t107+coef_J22.*t116+coef_J43.*t104.*2.0+coef_J23.*t125+coef_J44.*t113.*2.0+coef_J55.*t102+coef_J45.*t122.*2.0+coef_J56.*t111+coef_J57.*t120,coef_J18+coef_J21.*t108+coef_J22.*t117+coef_J23.*t126+coef_J50.*t103+coef_J51.*t112+coef_J52.*t121];
mt2 = [coef_J19+coef_J21.*t109+coef_J22.*t118+coef_J23.*t127+coef_J50.*t104+coef_J55.*t103+coef_J51.*t113+coef_J56.*t112+coef_J52.*t122+coef_J57.*t121,coef_J20+coef_J21.*t128+coef_J22.*t130+coef_J23.*t132+coef_J55.*t104+coef_J56.*t113+coef_J57.*t122,coef_J24+coef_J21.*t129+coef_J22.*t131+coef_J23.*t133,coef_J37.*4.0+coef_J43.*t105.*2.0+coef_J44.*t114.*2.0+coef_J45.*t123.*2.0,coef_J38.*3.0+coef_J43.*t106.*2.0+coef_J50.*t105+coef_J44.*t115.*2.0+coef_J51.*t114+coef_J45.*t124.*2.0+coef_J52.*t123,coef_J39.*3.0+coef_J43.*t107.*2.0+coef_J44.*t116.*2.0+coef_J55.*t105+coef_J45.*t125.*2.0+coef_J56.*t114+coef_J57.*t123,coef_J40.*2.0+coef_J43.*t108.*2.0+coef_J50.*t106+coef_J44.*t117.*2.0+coef_J51.*t115+coef_J45.*t126.*2.0+coef_J52.*t124,coef_J41.*2.0+coef_J43.*t109.*2.0+coef_J50.*t107+coef_J55.*t106+coef_J44.*t118.*2.0+coef_J51.*t116+coef_J56.*t115+coef_J45.*t127.*2.0+coef_J52.*t125+coef_J57.*t124];
mt3 = [coef_J42.*2.0+coef_J55.*t107+coef_J43.*t128.*2.0+coef_J56.*t116+coef_J44.*t130.*2.0+coef_J45.*t132.*2.0+coef_J57.*t125,coef_J46.*2.0+coef_J43.*t129.*2.0+coef_J44.*t131.*2.0+coef_J45.*t133.*2.0,coef_J47+coef_J50.*t108+coef_J51.*t117+coef_J52.*t126,coef_J48+coef_J50.*t109+coef_J55.*t108+coef_J51.*t118+coef_J56.*t117+coef_J52.*t127+coef_J57.*t126,coef_J49+coef_J55.*t109+coef_J56.*t118+coef_J50.*t128+coef_J51.*t130+coef_J52.*t132+coef_J57.*t127,coef_J53+coef_J50.*t129+coef_J51.*t131+coef_J52.*t133,coef_J54+coef_J55.*t128+coef_J56.*t130+coef_J57.*t132,coef_J58+coef_J55.*t129+coef_J56.*t131+coef_J57.*t133];
obj = [mt1,mt2,mt3];
end