/* ============================================================
   MedScope / SeqLab 本地预览服务
   用途：让页面从 http://localhost 加载，避免 file:// 下
        浏览器对多模态接口的跨域拦截。
   用法：node serve.cjs  （或双击「启动本地服务.cmd」）
   ============================================================ */
const http = require('node:http');
const fs = require('node:fs');
const path = require('node:path');

const ROOT = __dirname;
const PORT = Number(process.env.PORT) || 8000;

const MIME = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.png': 'image/png',
  '.jpg': 'image/jpeg',
  '.jpeg': 'image/jpeg',
  '.svg': 'image/svg+xml',
  '.ico': 'image/x-icon',
  '.md': 'text/markdown; charset=utf-8',
  '.fasta': 'text/plain; charset=utf-8',
  '.txt': 'text/plain; charset=utf-8',
};

const INDEX_HTML = `<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<title>bioinfo2026 · 本地预览</title>
<style>
  body{margin:0;min-height:100vh;display:grid;place-items:center;background:#0b1015;color:#c8d6e2;
       font-family:"Segoe UI","Microsoft YaHei",system-ui,sans-serif}
  .box{width:min(560px,calc(100vw - 48px));padding:32px 34px;background:#111a23;border:1px solid #243140;border-radius:16px;
       box-shadow:0 30px 80px -40px #000}
  h1{margin:0 0 6px;font-size:18px;letter-spacing:.04em}
  p.sub{margin:0 0 22px;font-size:12px;color:#5d6e7e}
  a{display:flex;align-items:center;justify-content:space-between;gap:12px;padding:14px 16px;margin-bottom:10px;
    border:1px solid #243140;border-radius:10px;text-decoration:none;color:#c8d6e2;background:#0d141a;transition:.15s}
  a:hover{border-color:#22d3ee66;background:#0f1a22;color:#e6f6fb}
  a b{font-weight:600;font-size:13.5px}
  a span{font-size:11px;color:#5d6e7e}
  .tip{margin-top:18px;font-size:11px;line-height:1.7;color:#5d6e7e}
</style></head><body>
<div class="box">
  <h1>期末项目 · 本地预览</h1>
  <p class="sub">http://localhost:${PORT}</p>
  <a href="/02-Imaging/index.html"><b>MedScope · 医学影像工作站</b><span>MPR / 3D / 影像组学 / AI 阅片</span></a>
  <a href="/01-SeqLab/index.html"><b>SeqLab · 序列分析工作台</b><span>ORF / 酶切图谱 / 密码子 / 比对</span></a>
  <div class="tip">关闭启动服务的命令行窗口即可停止服务。<br>通过 localhost 访问时，AI 面板调用多模态接口不会被浏览器跨域拦截。</div>
</div></body></html>`;

const server = http.createServer((req, res) => {
  let rel;
  try {
    rel = decodeURIComponent((req.url || '/').split('?')[0]);
  } catch (e) {
    res.writeHead(400); res.end('400'); return;
  }

  if (rel === '/' || rel === '/index.html') {
    res.writeHead(200, { 'Content-Type': 'text/html; charset=utf-8', 'Cache-Control': 'no-store' });
    res.end(INDEX_HTML);
    return;
  }

  const file = path.join(ROOT, rel);
  const rootWithSep = ROOT.endsWith(path.sep) ? ROOT : ROOT + path.sep;
  if (!file.startsWith(rootWithSep)) {
    res.writeHead(403, { 'Content-Type': 'text/plain; charset=utf-8' });
    res.end('403 Forbidden');
    return;
  }

  let target = file;
  try {
    if (fs.statSync(target).isDirectory()) target = path.join(target, 'index.html');
  } catch (e) {
    res.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' });
    res.end('404 Not Found: ' + rel);
    return;
  }

  if (!fs.existsSync(target)) {
    res.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' });
    res.end('404 Not Found: ' + rel);
    return;
  }

  const ext = path.extname(target).toLowerCase();
  res.writeHead(200, {
    'Content-Type': MIME[ext] || 'application/octet-stream',
    'Cache-Control': 'no-store',
  });
  fs.createReadStream(target).pipe(res);
});

server.listen(PORT, '127.0.0.1', () => {
  console.log('');
  console.log('  bioinfo2026 local preview');
  console.log('  ---------------------------------------------');
  console.log(`  Home           : http://localhost:${PORT}/`);
  console.log(`  Medical image  : http://localhost:${PORT}/02-Imaging/index.html`);
  console.log(`  Sequence lab   : http://localhost:${PORT}/01-SeqLab/index.html`);
  console.log('');
  console.log('  Press Ctrl+C or close this window to stop.');
  console.log('');
});
