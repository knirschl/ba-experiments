SUB=".generax_pick"

find . -name "*$SUB*" -exec rename "$SUB" "" {} \;

for file in ../../output/families/*/runs/F81/generax_eval_run/results/*; do
	if [[ "$file" == *"$SUB"* ]]; then
		rm -rf "${file%$SUB}"
		mv -v "$file" "${file%$SUB}"
	fi
done
