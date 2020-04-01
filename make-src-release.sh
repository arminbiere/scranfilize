#!/bin/sh
VERSION=`cat VERSION`
NAME=scranfilize-$VERSION
rm -rf /tmp/$NAME
mkdir /tmp/$NAME
cp -p \
configure \
LICENSE \
make-config \
makefile.in \
README.md \
scranfilize.c \
test.sh \
VERSION \
/tmp/$NAME/
sed \
  -i \
  -e "s,^GITID=.*,GITID=`git show|head -1|awk '{print $2}'`," \
/tmp/$NAME/make-config
mkdir /tmp/$NAME/cnfs
cp -p cnfs/*.cnf /tmp/$NAME/cnfs
cd /tmp/
rm -f $NAME.tar.xz
tar -cJf $NAME.tar.xz $NAME
ls -l /tmp/$NAME.tar.xz
