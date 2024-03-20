var1=$1
sed "s/.*rNstar.*/      rNstar = $var1/" input_orig.nlm > input.nlm
#awk 'NR==52 {$0=$mytext} 1' input.nlm
#awk '{ if (NR == 52) print "HI"; else print $0}' input.nlm input_orig.nlm
#cat input_orig.nlm | sed -e "s/the_original_line/the_new_line/" > temp_file.nlm
#sed '3/.*/HO/' input.nlm > new_file.nlm
