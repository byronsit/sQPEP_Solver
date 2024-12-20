function obj = coef_J_pure_pTop_sQPEP_merge_s(in1,in2,in3)
nv_1 = in3(1,:);
nv_2 = in3(2,:);
nv_3 = in3(3,:);
pa_1 = in1(1,:);
pa_2 = in1(2,:);
pa_3 = in1(3,:);
pb_1 = in2(1,:);
pb_2 = in2(2,:);
pb_3 = in2(3,:);
t2 = nv_1.*pa_1;
t3 = nv_2.*pa_2;
t4 = nv_3.*pa_3;
t5 = nv_1.*pb_1;
t6 = nv_2.*pb_2;
t7 = nv_3.*pb_3;
t8 = nv_1.*pa_2.*2.0;
t9 = nv_2.*pa_1.*2.0;
t10 = nv_1.*pa_3.*2.0;
t11 = nv_3.*pa_1.*2.0;
t12 = nv_2.*pa_3.*2.0;
t13 = nv_3.*pa_2.*2.0;
t14 = -t2;
t15 = -t9;
t16 = -t3;
t17 = -t11;
t18 = -t13;
t19 = -t4;
t20 = t8+t9;
t21 = t10+t11;
t22 = t12+t13;
t23 = t2+t3+t4;
t24 = t5+t6+t7;
t25 = t8+t15;
t26 = t10+t17;
t27 = t12+t18;
t28 = t2+t3+t19;
t29 = t2+t4+t16;
t30 = t3+t4+t14;
mt1 = [t23.^2,t23.*t27.*2.0,t23.*t26.*-2.0,t23.*t25.*2.0,t23.*t30.*-2.0+t27.^2,t20.*t23.*2.0-t26.*t27.*2.0,t21.*t23.*2.0+t25.*t27.*2.0,t23.*t29.*-2.0+t26.^2,t22.*t23.*2.0-t25.*t26.*2.0,t23.*t28.*-2.0+t25.^2,nv_1.*t23.*2.0,nv_2.*t23.*2.0,nv_3.*t23.*2.0,t23.*t24.*-2.0,t27.*t30.*-2.0,t20.*t27.*2.0+t26.*t30.*2.0,t21.*t27.*2.0-t25.*t30.*2.0,t20.*t26.*-2.0-t27.*t29.*2.0,t20.*t25.*2.0-t21.*t26.*2.0+t22.*t27.*2.0,t21.*t25.*2.0-t27.*t28.*2.0,nv_1.*t27.*2.0,nv_2.*t27.*2.0,nv_3.*t27.*2.0,t24.*t27.*-2.0,t26.*t29.*2.0,t22.*t26.*-2.0-t25.*t29.*2.0,t22.*t25.*2.0+t26.*t28.*2.0,nv_1.*t26.*-2.0,nv_2.*t26.*-2.0];
mt2 = [nv_3.*t26.*-2.0,t24.*t26.*2.0,t25.*t28.*-2.0,nv_1.*t25.*2.0,nv_2.*t25.*2.0,nv_3.*t25.*2.0,t24.*t25.*-2.0,t30.^2,t20.*t30.*-2.0,t21.*t30.*-2.0,t29.*t30.*2.0+t20.^2,t20.*t21.*2.0-t22.*t30.*2.0,t28.*t30.*2.0+t21.^2,nv_1.*t30.*-2.0,nv_2.*t30.*-2.0,nv_3.*t30.*-2.0,t24.*t30.*2.0,t20.*t29.*-2.0,t20.*t22.*2.0-t21.*t29.*2.0,t21.*t22.*2.0-t20.*t28.*2.0,nv_1.*t20.*2.0,nv_2.*t20.*2.0,nv_3.*t20.*2.0,t20.*t24.*-2.0,t21.*t28.*-2.0,nv_1.*t21.*2.0,nv_2.*t21.*2.0,nv_3.*t21.*2.0,t21.*t24.*-2.0,t29.^2,t22.*t29.*-2.0,t28.*t29.*2.0+t22.^2,nv_1.*t29.*-2.0,nv_2.*t29.*-2.0,nv_3.*t29.*-2.0,t24.*t29.*2.0,t22.*t28.*-2.0,nv_1.*t22.*2.0];
mt3 = [nv_2.*t22.*2.0,nv_3.*t22.*2.0,t22.*t24.*-2.0,t28.^2,nv_1.*t28.*-2.0,nv_2.*t28.*-2.0,nv_3.*t28.*-2.0,t24.*t28.*2.0,nv_1.^2,nv_1.*nv_2.*2.0,nv_1.*nv_3.*2.0,nv_1.*t24.*-2.0,nv_2.^2,nv_2.*nv_3.*2.0,nv_2.*t24.*-2.0,nv_3.^2,nv_3.*t24.*-2.0,t24.^2];
obj = [mt1,mt2,mt3];
end
