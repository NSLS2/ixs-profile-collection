c = get_config()
c.StoreMagics.autorestore = True
c.InteractiveShellApp.extensions = ['pyOlog.cli.ipy']
c.TerminalIPythonApp.log_datefmt = '%Y-%m-%d %H:%M:%S'
c.TerminalInteractiveShell.show_rewritten_input = True
c.TerminalInteractiveShell.autocall = 2
