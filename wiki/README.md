# Wiki source

These pages are the source of the
[pyfdstools wiki](https://github.com/johodges/pyfdstools/wiki). They live
here so that documentation changes are reviewed alongside the code
changes they describe.

## Publishing

The GitHub wiki is a separate git repository. To push these pages to it:

```bash
git clone https://github.com/johodges/pyfdstools.wiki.git /tmp/pyfdstools.wiki
cp wiki/*.md /tmp/pyfdstools.wiki/
rm /tmp/pyfdstools.wiki/README.md          # this file is not a wiki page
cd /tmp/pyfdstools.wiki
git add -A
git commit -m "Update wiki from the main repository"
git push
```

The file name becomes the page name, so `Quick-Start.md` is reachable at
`.../wiki/Quick-Start`. `Home.md` is the landing page and `_Sidebar.md`
is the navigation panel shown beside every page.

## Editing

Every complete code example on these pages was executed against the FDS
cases bundled in `pyfdstools/examples/data` before being committed. Run
anything you add before committing it — a documentation example that
does not run is worse than no example.

Internal links use the bare page name, `[Quick Start](Quick-Start)`, not
a path or a URL.
