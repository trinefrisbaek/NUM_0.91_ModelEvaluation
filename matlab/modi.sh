var1=$1
sed "s/.*rNstar.*/      rNstar = $var1/" ../input/input_orig.nlm > ../input/input.nlm
