incomings = input('trail # & frame # in an array\n');
if length(incomings) == 2
    wormspeed{incomings(1)}(1:incomings(2)) = -wormspeed{incomings(1)}(1:incomings(2));
elseif length(incomings) == 3
    wormspeed{incomings(1)}(incomings(2):incomings(3)) = -wormspeed{incomings(1)}(incomings(2):incomings(3));
end