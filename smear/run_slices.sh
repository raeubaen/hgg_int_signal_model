echo "sigma,errsigma,_,__"
source slices.sh $1 2> /dev/null | grep Sigma | grep -oP '(?<=Sigma).*' | tr -s ' ' | sed -E 's/^[[:space:]]+//' | sed -E 's\ \,\g'

