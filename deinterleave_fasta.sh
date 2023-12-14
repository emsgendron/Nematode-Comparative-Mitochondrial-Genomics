#!/bin/bash

sed -i -r 's/(>.*)/@\1@/' $1
wait
sed ':a;N;$!ba;s/\n//g' $1 >tmp1
wait
sed 's/@/\n/g' tmp1 >tmp2
sed '1d' tmp2 >tmp3
wait
mv tmp3 $1
rm tmp1
rm tmp2
