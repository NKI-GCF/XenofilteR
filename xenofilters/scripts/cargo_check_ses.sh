

cargo check --message-format short 2>&1 | grep error | python scripts/cargo_check_short_parser.py | tee /dev/stderr | xclip -selection clipboard

python scripts/find_conflicts.py --comments |  xclip -selection clipboard

git ls-files | grep -E '(\.rs|\.toml)' | ctags --fields=+n -L /dev/stdin
#git ls-files | grep -Ev '(\.(md|txt|sh)$|^todo/)' | ctags -L /dev/stdin  --fields=+n
#ctags --excmd=nr -R src

cat tags | python ctags_parser.py --skip-tests /dev/stdin |  xclip -selection clipboard

#cargo check --message-format short 2>&1 | sed -n -r 's/^.*`&*([^`]+)`.*$/^\1\t/; /[<>:]/!p' | sort | uniq  | grep -f /dev/stdin tags | xclip -selection clipboard
cargo check --message-format short 2>&1 | grep error | grep -Eo '`[^`]+`' | grep -v '[<>:]' | tr -d '`&' | sort | uniq > /tmp/check_matches.txt

sed 's/^/^/; s/$/[^a-zA-z0-9]/' /tmp/check_matches.txt  | grep -f /dev/stdin tags | xclip -selection clipboard
#You are a rust expert. Do not hallucinate. List per source file what more to use.

git diff master | python filter_diff.py --matches /tmp/check_matches.txt


python ctags_parser.py tags | grep -f <((echo [./];cargo check --message-format short 2>&1 | sed -n -r 's/^.*`&*([^`]+)`.*$/ \1$/; /[<>:]/!p' | sort | uniq)) | xclip -selection clipboard

cargo dupes -s --min-nodes 50 | python scripts/parse_fdupes_output.py --code --compact --prompt --no-comments | xclip -selection clipboard

python scripts/cargo_parser.py > /tmp/check_matches.txt

cat /tmp/check_matches.txt | python scripts/core_orient.py

cat /tmp/check_matches.txt | python scripts/core_orient.py --file-line | python scripts/fl_to_gd.py -U0

cargo check --message-format short 2>&1| grep error | xclip -selection clipboard


cargo check --message-format short 2>&1| grep error | python scripts/cargo_check_short_parser.py | xclip -selection clipboard

cat tags | python ctags_parser.py --skip-tests /dev/stdin |  xclip -selection clipboard

#cargo check --message-format short 2>&1 | sed -n -r 's/^.*`&*([^`]+)`.*$/^\1\t/; /[<>:]/!p' | sort | uniq  | grep -f /dev/stdin tags | xclip -selection clipboard
cargo check --message-format short 2>&1 | grep error | grep -Eo '`[^`]+`' | grep -v '[<>:]' | tr -d '`&' | sort | uniq > /tmp/check_matches.txt

sed 's/^/^/; s/$/[^a-zA-z0-9]/' /tmp/check_matches.txt  | grep -f /dev/stdin tags | xclip -selection clipboard
#You are a rust expert. Do not hallucinate. List per source file what more to use.

git diff master | python scripts/filter_diff.py --matches /tmp/check_matches.txt


python ctags_parser.py tags | grep -f <((echo [./];cargo check --message-format short 2>&1 | sed -n -r 's/^.*`&*([^`]+)`.*$/ \1$/; /[<>:]/!p' | sort | uniq)) | xclip -selection clipboard
