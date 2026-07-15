cargo check --message-format short 2>&1| grep error | xclip -selection clipboard

ctags --excmd=nr -R src

cat tags | python ctags_parser.py --skip-tests /dev/stdin |  xclip -selection clipboard

#cargo check --message-format short 2>&1 | sed -n -r 's/^.*`&*([^`]+)`.*$/^\1\t/; /[<>:]/!p' | sort | uniq  | grep -f /dev/stdin tags | xclip -selection clipboard
cargo check --message-format short 2>&1 | grep error | grep -Eo '`[^`]+`' | grep -v '[<>:]' | tr -d '`&' | sort | uniq > /tmp/check_matches.txt

sed 's/^/^/; s/$/[^a-zA-z0-9]/' /tmp/check_matches.txt  | grep -f /dev/stdin tags | xclip -selection clipboard
#You are a rust expert. Do not hallucinate. List per source file what more to use.

git diff master | python filter_diff.py --matches /tmp/check_matches.txt


python ctags_parser.py tags | grep -f <((echo [./];cargo check --message-format short 2>&1 | sed -n -r 's/^.*`&*([^`]+)`.*$/ \1$/; /[<>:]/!p' | sort | uniq)) | xclip -selection clipboard
