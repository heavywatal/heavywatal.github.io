module github.com/heavywatal/heavywatal.github.io

go 1.14

require (
	github.com/heavywatal/hugo-mod-common v0.0.0-20200815155143-dd3a158e3dfa // indirect
	github.com/heavywatal/hugo-theme-nonblog v0.0.0-20200815160457-faae40708a16
)

replace (
  github.com/heavywatal/hugo-mod-common => ./themes/hugo-mod-common
  github.com/heavywatal/hugo-theme-nonblog => ./themes/hugo-theme-nonblog
)
