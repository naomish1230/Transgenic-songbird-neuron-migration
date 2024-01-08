% RENUMBER CELL IDs SO THAT THEY GO IN ASCENDING ORDER (useful if combined multiple birds or deleted cells)
function  [current] = RenumberCells(current)

newN = [];
count = 1;
for n = 1:height(current)-1
if current.N(n) == current.N(n+1)
newN(n) = count;
newN(n+1) = count;
elseif current.N(n) ~= current.N(n+1)
newN(n) = count;
newN(n+1) = count+1;
count = count+1;
end
end
newN = newN';
current.N = newN;

end
