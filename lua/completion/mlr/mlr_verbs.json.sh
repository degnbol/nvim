#!/usr/bin/env zsh
# Generate mlr_verbs.json from the zsh completion data.
# Reads from: ~/dotfiles/miller/completion/

set -euo pipefail
cd ${0:A:h}
local MILLER_DIR=~/dotfiles/miller/completion

if [[ ! -d "$MILLER_DIR/verb" ]]; then
    echo "Error: $MILLER_DIR/verb not found. Run RUNME.zsh in miller/completion first." >&2
    exit 1
fi

# Build JSON with jq
local verbs_json='{}'
while IFS=: read -r name desc; do
    desc="${desc## }"
    (( ${#desc} > 100 )) && desc="${desc:0:97}..."
    verbs_json=$(printf '%s' "$verbs_json" | jq --arg k "$name" --arg v "$desc" '. + {($k): $v}')
done < "$MILLER_DIR/descs.help"

local field_flags_json='{}'
local verb_flags_json='{}'

for opt_file in "$MILLER_DIR"/verb/*.help.opt(N); do
    local verb="${${opt_file:t}%.help.opt}"
    local -a ff=()
    local verb_flags_arr='[]'

    while IFS= read -r comp; do
        if [[ "$comp" =~ '^(-[-a-zA-Z]+)(.*)' ]]; then
            local flag="${match[1]}"
            local rest="${match[2]}"
            local desc=""
            if [[ "$rest" == '['* ]]; then
                desc="${rest#\[}"
                desc="${desc%\]}"
                desc="${desc#\{left-file-fields\} }"
                (( ${#desc} > 80 )) && desc="${desc:0:77}..."
            fi
            verb_flags_arr=$(printf '%s' "$verb_flags_arr" | jq --arg f "$flag" --arg d "$desc" '. + [[$f, $d]]')

            # Detect field-taking flags (same patterns as _mlr.sh)
            if [[ "$comp" =~ '^-[-a-zA-Z]+\[(\{[a-z],[a-z],[a-z]\}|\{comma-separated|\{one or more comma-separated|Field name)' ]]; then
                ff+=("$flag")
            elif [[ "$comp" =~ '^-[-a-zA-Z]+\[\{left-file-fields\}' ]]; then
                ff+=("$flag")
            fi
        fi
    done < "$opt_file"

    verb_flags_json=$(printf '%s' "$verb_flags_json" | jq --arg k "$verb" --argjson v "$verb_flags_arr" '. + {($k): $v}')

    if (( ${#ff} > 0 )); then
        field_flags_json=$(printf '%s' "$field_flags_json" | jq \
            --arg k "$verb" \
            --argjson v "$(printf '%s\n' "${ff[@]}" | jq -R . | jq -s .)" \
            '. + {($k): $v}')
    fi
done

jq -n \
    --argjson verbs "$verbs_json" \
    --argjson field_flags "$field_flags_json" \
    --argjson verb_flags "$verb_flags_json" \
    '{
        verbs: $verbs,
        field_flags: $field_flags,
        verb_flags: $verb_flags,
        positional_field_verbs: ["group-by", "label", "rename", "sec2gmt", "sec2gmtdate"]
    }' > mlr_verbs.json
