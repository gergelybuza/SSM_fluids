function LnB = get_LnB(obj,LAM,k,om_r)
% L minus i*om*B
LnB = obj.get_L(LAM,k) - 1i*om_r*obj.B;
end