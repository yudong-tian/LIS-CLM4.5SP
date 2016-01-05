
file=/discover/nobackup/projects/nca/ytian/CESM_INPUT_DATA/atm/datm7/domain.T62.050609.nc
ncdump -h $file | egrep '(int |float |double ).*\(.*\)' > tmp.vars 

while read var; do 

  cmd=`echo $var | sed 's/[,();]/ /g'`
  varname=`echo $cmd | awk '{print \$2}'` 
  #echo $varname 

  nargs=`echo $cmd |wc -w`
  let ndims=nargs-2

  cat > ${varname}.ctl <<EOF
DSET ^$varname.1gd4r 
undef 1.e+36
options big_endian 
EOF

dim_sizes=`./read_netcdf $file $cmd` 
xdim=`echo $dim_sizes |cut -d' ' -f1`
ydim=`echo $dim_sizes |cut -d' ' -f2`
zdim=`echo $dim_sizes |cut -d' ' -f3`
tdim=`echo $dim_sizes |cut -d' ' -f4`

zlev=2     # if multiple zlevs, display the 2nd level, as 1st level may be zero
if [ -z "$zdim" ]; then 
  zdim=1
  zlev=1
fi
if [ -z "$tdim" ]; then 
  tdim=1
fi

gvarname=`echo $varname |cut -c1-15`

  cat >> ${varname}.ctl <<EOF1
XDEF $xdim   LINEAR    1   1
YDEF $ydim   LINEAR    1   1
ZDEF $zdim   LINEAR    1 1 
TDEF $tdim   LINEAR 00Z01Jan1901 1mo
VARS 1
$gvarname $zdim 99 $varname
ENDVARS
EOF1

grads -bl <<EOF3
open ${varname}.ctl 
set gxout shaded
set grads off
set z $zlev
d ${gvarname} 
draw xlab ni
draw ylab nj
cbar 
draw title $varname
gxyat -x 1000 -y 800 ${varname}.png
quit
EOF3

done < tmp.vars

# generate html page 

cat > index.html <<EOF4
<html>
<head></head>
<body> 
<table> 
<tr>
EOF4

ncols=4
i=0

sort -f  -k 2 tmp.vars |
{ while read var; do
  ((i++))
  echo i=$i
  cmd=`echo $var | sed 's/[,();]/ /g'`
  varname=`echo $cmd | awk '{print \$2}'`
  cat >> index.html <<EOFx
<td><a href="$varname.png"><img src="$varname.png" width=300></a></td>
EOFx

 if (( $i % $ncols == 0 )); then
   i=0
   cat >> index.html <<EOFy
</tr><tr>
EOFy
 fi

done
}

cat >> index.html <<EOF5
</tr>
</table>
</body>
</html>
EOF5

