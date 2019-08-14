#!/usr/bin/env bash
if [ "x$BASH" = x ] || [ ! "$BASH_VERSINFO" ] || [ "$BASH_VERSINFO" -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue
unset CDPATH

TestDir=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
DunovoDir=$(dirname "$TestDir")
BfxDir="$(dirname "$DunovoDir")/nick-bfx"
RefName="overlap.ref.fa"
Fq1Name="overlap.reads_1.fq"
Fq2Name="overlap.reads_2.fq"
FamiliesName="overlap.families.tsv"
MsaName="overlap.families.msa.tsv"
Sscs1Name="overlap.sscs_1.fa"
Sscs2Name="overlap.sscs_2.fa"
BamName="overlap.sscs.bam"
TmpOutputs="$FamiliesName $Sscs1Name $Sscs2Name"
Outputs="$TmpOutputs $RefName $Fq1Name $Fq2Name $MsaName $BamName"

Usage="Usage: \$ $(basename "$0") -d
       \$ $(basename "$0") [-k] overlap.align.txt
Options:
-d: Delete existing output files.
-k: Keep intermediate temporary files."

function main {

  # Get arguments.
  keep=
  delete=
  while getopts "dkh" opt; do
    case "$opt" in
      d) delete="true";;
      k) keep="true";;
      [h?]) fail "$Usage";;
    esac
  done
  align="${@:$OPTIND:1}"

  if ! ([[ "$align" ]] || [[ "$delete" ]]); then
    fail "$Usage"
  fi

  if [[ "$delete" ]]; then
    delete_outputs "$Outputs"
    return "$?"
  fi

  if existing_outputs "$Outputs"; then
    return 1
  fi

  missing_prereq=
  for script in "$TestDir/parse-test-align.py" "$DunovoDir/make-families.sh" \
                "$DunovoDir/align-families.py" "$BfxDir/align-bowtie.sh"; do
    if ! [[ -x "$script" ]]; then
      missing_prereq="true"
      echo "Error: Could not find or execute script $script" >&2
    fi
  done
  if [[ "$missing_prereq" ]]; then
    return 1
  fi

  barlen=$("$TestDir/parse-test-align.py" --print-barlen --ref "$TestDir/$RefName" \
           --fq1 "$TestDir/$Fq1Name" --fq2 "$TestDir/$Fq2Name" "$align")

  "$DunovoDir/make-families.sh" -t "$((barlen/2))" "$TestDir/$Fq1Name" "$TestDir/$Fq2Name" \
    > "$TestDir/$FamiliesName"

  "$DunovoDir/align-families.py" "$TestDir/$FamiliesName" > "$TestDir/$MsaName"

  "$DunovoDir/make-consensi.py" "$TestDir/$MsaName" --min-reads 3 --qual 20 \
    --sscs1 "$TestDir/$Sscs1Name" --sscs2 "$TestDir/$Sscs2Name"

  "$BfxDir/align-bowtie.sh" -c -b "$TestDir/$BamName" "$TestDir/$RefName" \
    "$TestDir/$Sscs1Name" "$TestDir/$Sscs2Name"

  if ! [[ "$keep" ]]; then
    for filename in $TmpOutputs; do
      rm "$TestDir/$filename"
    done
  fi
}


function delete_outputs {
  outputs="$1"
  existing=
  echo "About to delete these files in directory $TestDir:"
  for filename in $outputs; do
    if [[ -e "$TestDir/$filename" ]]; then
      echo "  $filename"
      existing="true"
    fi
  done
  if ! [[ "$existing" ]]; then
    echo "No existing files found."
    return 0
  fi
  echo -n "Type \"yes\" to proceed: "
  read response
  if [[ "$response" != "yes" ]]; then
    echo "Aborting.."
    return 1
  fi
  for filename in $outputs; do
    path="$TestDir/$filename"
    if [[ -f "$path" ]]; then
      rm "$path"
    fi
  done
}


function existing_outputs {
  outputs="$1"
  existing=
  for filename in $outputs; do
    path="$TestDir/$filename"
    if [[ -e "$path" ]]; then
      existing="$existing $filename"
    fi
  done
  if [[ "$existing" ]]; then
    echo "Error: The following files already exist in $TestDir:" >&2
    for filename in $existing; do
      echo "  $filename"
    done
    return 0
  else
    return 1
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
