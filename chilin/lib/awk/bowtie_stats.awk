BEGIN{total=0; map=0; a=0;b=0}
{
if (/^[^@]/){
total+=1 
if ($2!="4") {
    map+=1 
    ur[$1] += 1
    ul[$2":"$3":"$4] += 1
}
}
}
END{
for (urr in ur)
a+=1
for (ull in ul)
b+=1
print total
print map
print a
print b
}
