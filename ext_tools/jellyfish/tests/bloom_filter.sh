#! /bin/sh

mkdir -p tests-data; cd tests-data
. ../compat.sh

# sort -k2,2 > ${pref}.md5sum <<EOF
# EOF

$JF count --bf-size 10M --bf-fp 0.001 -t $nCPUs -o ${pref}_10m.jf -s 1M -m 40 seq10m.fa
$JF histo ${pref}_10m.jf > ${pref}_10m.histo
if bc < /dev/null 2> /dev/null; then
    COLLISIONS=$(cut -d\  -f2 ${pref}_10m.histo | paste -sd+ - | bc)
    [ $((COLLISIONS > 10000)) = 0 ] || {
        echo >&2 "Too many collisions"
        false
    }
else
    echo > "bc missing: Skip collisions test"
fi

$JF count --bf-size 3M --bf-fp 0.001 -t $nCPUs -o ${pref}_3m.jf -s 1M -m 40 seq1m_0.fa seq1m_1.fa seq1m_0.fa seq1m_2.fa
$JF histo ${pref}_3m.jf > ${pref}_3m.histo

if bc < /dev/null 2> /dev/null; then
    COLLISIONS=$(cut -d\  -f2 ${pref}_3m.histo | paste -sd+ - | bc)
    [ $((COLLISIONS - 1000000 > 20000)) = 0 ] || {
        echo >&2 "Too many collisions"
        false
    }
else
    echo "bc missing: Skip collisions test"
fi
