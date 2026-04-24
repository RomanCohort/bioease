# 运行download_TCGA_data.R文件下载TCGA数据
# 再运行analysis.R文件进行分析
# 运行front.py文件进行可视化

import os
import sys
import subprocess
import shutil
import streamlit as st
from streamlit.runtime.scriptrunner import get_script_run_ctx
from streamlit.web import bootstrap


def ensure_streamlit_context():
    """Fallback for accidental direct execution (including stale packaged entrypoints)."""
    if get_script_run_ctx() is not None:
        return
    if os.environ.get('BIOEASE_FRONT_BOOTSTRAPPED') == '1':
        return
    os.environ['BIOEASE_FRONT_BOOTSTRAPPED'] = '1'
    bootstrap.run(os.path.abspath(__file__), False, [], {})
    raise SystemExit(0)


ensure_streamlit_context()

def get_app_root():
    if getattr(sys, 'frozen', False):
        return getattr(sys, '_MEIPASS', os.path.dirname(sys.executable))
    return os.path.dirname(os.path.abspath(__file__))


def ensure_results_dir(app_root):
    candidates = [
        os.path.join(os.getcwd(), 'results'),
        os.path.join(app_root, 'results'),
        os.path.join(os.environ.get('LOCALAPPDATA', os.path.expanduser('~')), 'bioease', 'results'),
    ]
    for path in candidates:
        try:
            os.makedirs(path, exist_ok=True)
            return path
        except PermissionError:
            continue
    raise PermissionError('无法创建 results 目录，请检查当前目录和用户目录权限。')


APP_ROOT = get_app_root()
RESULTS_DIR = ensure_results_dir(APP_ROOT)

st.title("BIOEASE(*'ω'*)")

st.markdown("""
这是一个生物信息学分析工具，用于下载TCGA数据并进行基因差异分析和生存分析。
**注意：** 本工具为练习项目，结果仅供参考。请确保R环境已安装并配置好相关包。
""")
st.caption(f"结果目录: {RESULTS_DIR}")

cancer = st.text_input("请输入你想研究的癌症类型", placeholder="例如：BRCA")
st.write("不知道癌症类型对应的代号？请看[这里](https://zhuanlan.zhihu.com/p/467229434)")
agree = st.checkbox("我确认本工具为练习项目，结果仅供参考")

# 自动检测 Rscript 可执行文件
def detect_rscript():
    # 优先使用系统路径
    exe = shutil.which('Rscript') or shutil.which('Rscript.exe')
    if exe:
        return exe
    # 在常见安装目录中查找
    base = r"C:\Program Files\R"
    if os.path.isdir(base):
        for name in os.listdir(base):
            candidate = os.path.join(base, name, 'bin', 'Rscript.exe')
            if os.path.isfile(candidate):
                return candidate
    # 还可以检查 Program Files (x86)
    base2 = r"C:\Program Files (x86)\R"
    if os.path.isdir(base2):
        for name in os.listdir(base2):
            candidate = os.path.join(base2, name, 'bin', 'Rscript.exe')
            if os.path.isfile(candidate):
                return candidate
    return None

detected_rscript = detect_rscript()
provided_rscript = st.text_input('Rscript 可执行路径（可选，若自动检测失败请在此粘贴完整路径）', value=detected_rscript or '')
rscript_exe = provided_rscript or detected_rscript or ''

if not rscript_exe:
    st.warning('未检测到 Rscript。请确保已安装 R 并将其 bin 目录加入 PATH，或在上方输入 Rscript 可执行文件完整路径。')
else:
    st.info(f'使用 Rscript: {rscript_exe}')

# 验证 Rscript 可执行性（避免 PermissionError）
def verify_rscript(path):
    """尝试运行 `Rscript --version` 来验证可执行权限。
    返回 (ok: bool, msg: str)。
    """
    if not path:
        return False, 'Rscript 路径为空'
    if not os.path.isfile(path):
        return False, f'找不到文件: {path}'
    try:
        env = os.environ.copy()
        env['PATH'] = env.get('PATH', '') + ';' + os.path.dirname(path)
        p = subprocess.run([path, '--version'], capture_output=True, text=True, check=True, env=env)
        return True, p.stdout or p.stderr or 'Rscript 可执行，未返回文本'
    except subprocess.CalledProcessError as e:
        # Rscript 返回非0，但说明可执行
        return True, e.stdout or e.stderr or 'Rscript 返回错误，但可执行'
    except PermissionError:
        return False, 'PermissionError: 无法执行该文件，访问被拒绝'
    except FileNotFoundError:
        return False, 'FileNotFoundError: 未找到可执行文件'
    except Exception as e:
        return False, f'其他错误: {str(e)}'

col1, col2, col3 = st.columns(3)

with col1:
    download_button = st.button("点击下载数据")
with col2:
    analysis_button = st.button("点击进行基因差异分析")
with col3:
    survival_button = st.button("点击进行生存分析")

if download_button:
    if not agree:
        st.warning("请先勾选确认框再生成。")
    elif not cancer.strip():
        st.error("请输入有效的癌症类型。")
    else:
        cancer_full = "TCGA-" + cancer.upper().strip()
        with st.spinner("正在下载数据...，下载时间可能较长"):
            ok, msg = verify_rscript(rscript_exe)
            if not ok:
                st.error(f"Rscript 验证失败：{msg}")
                st.markdown("请按下列步骤排查：")
                st.markdown("- 确认在页面顶部填写了 `Rscript.exe` 的完整路径；")
                st.markdown("- 在 PowerShell 中测试：")
                st.code(r"& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' --version")
                st.code(r"Get-Acl 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' | Format-List")
                st.code(r"icacls 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' /grant %USERNAME%:RX")
            else:
                try:
                    env = os.environ.copy()
                    env['PATH'] = env.get('PATH', '') + ';' + os.path.dirname(rscript_exe)
                    result = subprocess.run(
                        [rscript_exe, os.path.join(APP_ROOT, 'download_TCGA_data.R'), '--cancer', cancer_full],
                        capture_output=True,
                        text=True,
                        check=True,
                        env=env,
                        cwd=APP_ROOT,
                    )
                    st.success("数据下载完成！保存在当前目录下")
                    st.text_area("下载日志", result.stdout, height=100)
                except subprocess.CalledProcessError as e:
                    st.error(f"下载失败：{e.stderr or '无错误输出'}")
                    if e.stdout:
                        st.text_area("标准输出", e.stdout, height=100)
                except PermissionError:
                    st.error("PermissionError: 无法执行 Rscript（访问被拒绝）。请以管理员身份运行 Streamlit 或检查文件权限。")
                    st.code(r"Start-Process -Verb RunAs powershell;  # 在管理员权限下重启 PowerShell 并运行 streamlit run front.py")
                except FileNotFoundError:
                    st.error("Rscript未找到，请确保R已安装并在PATH中，或在页面上方填写可执行文件路径。")

if analysis_button:
    if not agree:
        st.warning("请先确认数据已下载完毕。")
    else:
        with st.spinner("正在进行基因差异分析..."):
            ok, msg = verify_rscript(rscript_exe)
            if not ok:
                st.error(f"Rscript 验证失败：{msg}")
                st.markdown("请首先按页面上方提示确认 Rscript 可执行性，然后重试。")
            else:
                try:
                    env = os.environ.copy()
                    env['PATH'] = env.get('PATH', '') + ';' + os.path.dirname(rscript_exe)
                    result = subprocess.run(
                        [rscript_exe, os.path.join(APP_ROOT, 'analysis.R')],
                        capture_output=True,
                        text=True,
                        check=True,
                        env=env,
                        cwd=APP_ROOT,
                    )
                    st.success("分析完成！火山图与基因计数矩阵保存在当前目录下")
                    st.text_area("分析日志", result.stdout, height=100)
                except subprocess.CalledProcessError as e:
                    st.error(f"分析失败：{e.stderr or '无错误输出'}")
                    if e.stdout:
                        st.text_area("标准输出", e.stdout, height=100)
                except PermissionError:
                    st.error("PermissionError: 无法执行 Rscript（访问被拒绝）。请以管理员身份运行 Streamlit 或检查文件权限。")
                except FileNotFoundError:
                    st.error("Rscript未找到，请确保R已安装并在PATH中，或在页面上方填写可执行文件路径。")

if survival_button:
    if not agree:
        st.warning("请先确认计数矩阵已生成。")
    else:
        with st.spinner("正在进行生存分析..."):
            ok, msg = verify_rscript(rscript_exe)
            if not ok:
                st.error(f"Rscript 验证失败：{msg}")
                st.markdown("请首先按页面上方提示确认 Rscript 可执行性，然后重试。")
            else:
                try:
                    env = os.environ.copy()
                    env['PATH'] = env.get('PATH', '') + ';' + os.path.dirname(rscript_exe)
                    result = subprocess.run(
                        [rscript_exe, os.path.join(APP_ROOT, 'upregulated_TAA_survival_analysis.R')],
                        capture_output=True,
                        text=True,
                        check=True,
                        env=env,
                        cwd=APP_ROOT,
                    )
                    st.success("生存分析完成！生存分析图保存在当前目录下")
                    st.text_area("生存分析日志", result.stdout, height=100)
                except subprocess.CalledProcessError as e:
                    st.error(f"生存分析失败：{e.stderr or '无错误输出'}")
                    if e.stdout:
                        st.text_area("标准输出", e.stdout, height=100)
                except PermissionError:
                    st.error("PermissionError: 无法执行 Rscript（访问被拒绝）。请以管理员身份运行 Streamlit 或检查文件权限。")
                except FileNotFoundError:
                    st.error("Rscript未找到，请确保R已安装并在PATH中，或在页面上方填写可执行文件路径。")