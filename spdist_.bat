echo d|xcopy %1\.Data %2\.Data /e /s
echo d|xcopy %1\.Prefs %2\.Prefs /e /s
echo d|xcopy %1\doc %2\doc /e /s
copy %1\S.dll %2\S.dll
copy %1\resample.chm %2\resample.chm
copy %1\DESCRIPTION %2\DESCRIPTION
copy %1\resample.pdf %2\resample.pdf
