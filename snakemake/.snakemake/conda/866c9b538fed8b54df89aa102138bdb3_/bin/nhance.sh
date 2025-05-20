#!/bin/bash

# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#            National Center for Biotechnology Information (NCBI)
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government do not place any restriction on its use or reproduction.
#  We would, however, appreciate having the NCBI and the author cited in
#  any work or product based on this material.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
# ===========================================================================
#
# File Name:  nhance.sh
#
# Author:  Jonathan Kans
#
# Version Creation Date:   06/15/20
#
# ==========================================================================

# external nquire shortcut extensions

xtra=""
ptrn=""
encd=""

if [ $# -gt 0 ]
then
  xtra="$1"
  case "$xtra" in
    -pathway | -gene-to-pathway | -litvar | -citmatch )
      shift
      if [ $# -gt 0 ]
      then
        ptrn="$1"
        encd=$( Escape "$ptrn" )
        shift
      else
        echo "ERROR: Missing $xtra argument" >&2
        exit 1
      fi
      ;;
    * )
      ;;
  esac
fi

# convert encoded entities, but convert &quot;, &apos, and &#39; (apostrophe) to space
# adapted from:
#   https://stackoverflow.com/questions/5929492/bash-script-to-convert-from-html-entities-to-characters
HTMLtoText () {
  LineOut=$1  # Parm 1= Input line
  # Replace external command: Line=$(sed 's/&amp;/\&/g; s/&lt;/\</g; 
  # s/&gt;/\>/g; s/&quot;/\"/g; s/&#39;/\'"'"'/g; s/&ldquo;/\"/g; 
  # s/&rdquo;/\"/g;' <<< "$Line") -- With faster builtin commands.
  LineOut="${LineOut//&nbsp;/ }"
  LineOut="${LineOut//&amp;/&}"
  LineOut="${LineOut//&lt;/<}"
  LineOut="${LineOut//&gt;/>}"
  LineOut="${LineOut//&quot;/ }"
  LineOut="${LineOut//&apos;/ }"
  LineOut="${LineOut//&#39;/ }"
  echo "${LineOut}"
}

if [ -n "$encd" ]
then
  case "$xtra" in
    -pathway )
      # Reactome:R-HSA-70171
      nquire -pugview data pathway "${encd}" XML
      ;;
    -gene-to-pathway )
      # 1956
      fst='{"download":"*","collection":"pathway","where":{"ands":[{"geneid":"'
      scd='"}]},"order":["taxname,asc"],"start":1,"limit":10000000,"downloadfilename":"GeneID_'
      thd='_pathway"}'
      jsn=$( echo "${fst}${encd}${scd}${encd}${thd}" )
      nquire -pubchem sdq/sdqagent.cgi -infmt json -outfmt xml -query "${jsn}"
      ;;
    -litvar )
      # rs11549407
      nquire -get "https://www.ncbi.nlm.nih.gov/research/bionlp/litvar" \
        "api/v1/entity/litvar/${encd}%23%23" |
      transmute -j2x
      ;;
    -citmatch )
      # clean encoded characters
      ptrn=$( HTMLtoText "$ptrn" )
      # echo "${ptrn}" >&2
      # "nucleotide sequences required for tn3 transposition immunity"
      nquire -get https://pubmed.ncbi.nlm.nih.gov/api/citmatch \
        -method heuristic -raw-text "$ptrn" |
      transmute -j2x |
      xtract -pattern opt \
        -if success -equals true -and result/count -eq 1 \
        -sep "\n" -element uids/pubmed
      ;;
    * )
      ;;
  esac

  exit 0
fi

