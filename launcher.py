import os
import socket
import sys
import threading
import time
import webbrowser
from datetime import datetime
from streamlit import config as st_config
from streamlit.web import bootstrap


def get_resource_root():
    if getattr(sys, 'frozen', False):
        return getattr(sys, '_MEIPASS', os.path.dirname(sys.executable))
    return os.path.dirname(os.path.abspath(__file__))


def get_log_file():
    base = os.environ.get('LOCALAPPDATA', os.path.expanduser('~'))
    log_dir = os.path.join(base, 'bioease', 'logs')
    os.makedirs(log_dir, exist_ok=True)
    return os.path.join(log_dir, 'launcher.log')


def log_line(message):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(get_log_file(), 'a', encoding='utf-8') as f:
        f.write(f'[{ts}] {message}\n')


def can_bind_local_port(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            sock.bind(("127.0.0.1", port))
            return True
        except OSError:
            return False


def select_port():
    # Prefer Streamlit default ports to avoid accidental 3000 dev-server behavior.
    candidates = [8501, 8502, 8503, 8504, 8505]
    for port in candidates:
        if can_bind_local_port(port):
            return port
    return 8501


def open_browser_later(url, delay_sec=1.2):
    def _open():
        time.sleep(delay_sec)
        webbrowser.open(url)

    threading.Thread(target=_open, daemon=True).start()


def main():
    script_path = os.path.join(get_resource_root(), 'front.py')
    if not os.path.isfile(script_path):
        raise FileNotFoundError(f'Cannot find Streamlit script: {script_path}')

    # Start Streamlit programmatically to avoid CLI target parsing against .exe argv.
    port = select_port()
    url = f'http://localhost:{port}'
    open_browser_later(url)
    log_line(f'script={script_path}')
    log_line(f'selected_port={port}')
    os.environ['STREAMLIT_GLOBAL_DEVELOPMENT_MODE'] = 'false'
    os.environ['STREAMLIT_SERVER_PORT'] = str(port)
    os.environ['STREAMLIT_SERVER_HEADLESS'] = 'true'
    os.environ['STREAMLIT_BROWSER_SERVER_PORT'] = str(port)
    flag_options = {
        'global.developmentMode': False,
        'server.port': port,
        'server.address': '127.0.0.1',
        'server.headless': True,
        'browser.serverPort': port,
        'browser.gatherUsageStats': False,
    }
    bootstrap.load_config_options(flag_options)
    # Force-set important options in case user/global config overrides flags.
    st_config.set_option('global.developmentMode', False)
    st_config.set_option('server.port', port)
    st_config.set_option('server.address', '127.0.0.1')
    st_config.set_option('server.headless', True)
    st_config.set_option('browser.serverPort', port)

    log_line(f"effective.global.developmentMode={st_config.get_option('global.developmentMode')}")
    log_line(f"effective.server.port={st_config.get_option('server.port')}")
    log_line(f"effective.browser.serverPort={st_config.get_option('browser.serverPort')}")
    log_line(f"effective.server.headless={st_config.get_option('server.headless')}")

    bootstrap.run(script_path, False, [], flag_options)


if __name__ == '__main__':
    main()
