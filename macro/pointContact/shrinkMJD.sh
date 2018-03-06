#/bin/bash
sed -i -e '1,36d' ev.dat
echo "">ev.new 
while read w1 w2 w3 w4 w5 w6 ;
do
  echo "$w1 $w2 $w3">>ev.new
  if  [ "$w1" != "0.00" ] 
  then
    echo "-$w1 $w2 $w3">>ev.new
  fi
done< ev.dat
