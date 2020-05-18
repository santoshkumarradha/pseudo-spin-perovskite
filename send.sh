#!/bin/bash  
git pull
git add .
git commit -m "$*"
git push --set-upstream origin master
