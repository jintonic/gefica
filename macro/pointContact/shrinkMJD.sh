#/bin/bash
sed -i -e '1,36d' "$1"
echo "">ev.new
while read w1 w2 w3 w4 w5 w6 ;
do
  echo "$w1 $w2 $w3">>ev.new
done< "$1"
