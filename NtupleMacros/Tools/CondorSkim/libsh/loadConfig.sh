#! /bin/bash

cfg=$1

var=
value=

#seems to cut off the last line for some reason, be sure to finish config file with a new line
while read -r line; do
	#remove comments, deliniated with a "#"
	validline="${line%%#*}"
	#split at the "=" sign and assign the first half to $var
	var="${validline%%=*}"
	#trim leading and trailing spaces from var
	while : #infinite loop to catch all spaces
	do
		case $var in
			' '*) var=${var#?};; #remove leading space
			*' ') var=${var%?};; #remove all trailing spaces
			*) break;; #once no more spaces, quit infite loop
		esac
	done
	#split at the "=" sign and assign the second half to $value
	value="${validline##*=}"
	#trim leading and trailing spaces from var
	while : #infinite loop to catch all spaces
	do
		case $value in
			' '*) value=${value#?};; #remove leading space
			*' ') value=${value%?};; #remove all trailing spaces
			*) break;; #once no more spaces, quit infite loop
		esac
	done

	#make sure that var exists
	if [ ${#var} -lt 1 ]; then
		#echo "var is too short"
		continue
	fi
	#make sure that var is a valid variable name || "$var" =~ *[!a-zA-Z0-9]* 
	if [[ "$var" != [a-zA-Z_]* || "$var" != *[a-zA-Z0-9_]* ]]; then
		#echo "not a valid variable name"
		continue
	fi
	#make sure that val exists
	if [ "${#value}" -lt 1 ]; then
		#echo "value is too short"
		continue
	fi

	#echo "$line:$validline:$var:$value"
	
	#finally save the variable
	eval "$var='$value'"
done < $cfg

