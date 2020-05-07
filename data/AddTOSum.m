function summa = AddToSum(summa, add)
if length(summa) < length(add)
        summa = [summa; zeros(length(add)-length(summa), 1)];
    elseif length(summa) > length(add)
        add = [add; zeros(-length(add)+length(summa), 1)];
    end
    summa = summa + add; 