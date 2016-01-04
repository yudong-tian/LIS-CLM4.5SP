
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

