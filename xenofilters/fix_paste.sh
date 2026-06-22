file="$1"

sed -i -r -n '
1b;
$b;
s/&gt;/\>/g;
s/&lt;/\</g;
s/&amp;/\&/g;
s/&quot;/\"/g;
s/&#39;/\x27/g;
s/&nbsp;/\xA0/g;
s/&Tab;/\t/g;
s/&NewLine;/\n/g;
s~^(</span>){1,2}<span>~~;
/^<span/!p;
' "$file"
