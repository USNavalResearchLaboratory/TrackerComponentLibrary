function l = get_scs_cone_dims(K)
l = K.f + K.l;
for i=1:length(K.q)
    l = l + K.q(i);
end
for i=1:length(K.s)
    l = l + get_sd_cone_size(K.s(i));
end
l = l + K.ep * 3;
l = l + K.ed * 3;
l = l + length(K.p) * 3;
end