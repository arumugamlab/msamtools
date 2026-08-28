# Contributing to msamtools

Contributions to msamtools are welcome.

## Reporting bugs and requesting features

Please use the GitHub issue tracker to report bugs, request features, or discuss changes in behavior.

When reporting a bug, please include:

- the msamtools version or commit used,
- the command that was run,
- the observed behavior,
- the expected behavior, and
- a minimal reproducible example where possible.

## Contributing code

Please make changes in a separate branch and submit them through a pull request.

Before opening a pull request, run:

```bash
autoreconf -fi
./configure
make
```

If you change command-line help text, regenerate the corresponding documentation before running the test suite:

```bash
scripts/update_cli_docs.sh
```

Then run:

```bash
make check
```

If documentation is modified, also run:

```bash
mkdocs build --strict
```

Please keep changes focused and follow the style of the surrounding code.

## Support

Questions about msamtools usage can be raised through the GitHub issue tracker.
