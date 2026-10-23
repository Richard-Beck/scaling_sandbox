import {spawn} from 'node:child_process';
import {mkdtempSync,writeFileSync,mkdirSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join,resolve,dirname} from 'node:path';
import {fileURLToPath} from 'node:url';
const methodFolder=dirname(fileURLToPath(import.meta.url));
const input=process.argv[2]||resolve(methodFolder,'../viewer/index.html');
const output=process.argv[3]||resolve(methodFolder,'browser_checks');
mkdirSync(output,{recursive:true});
const chrome=spawn(process.env.CHROME_BIN||'chromium',['--headless','--no-sandbox','--disable-dev-shm-usage','--disable-gpu','--remote-debugging-port=0','--remote-allow-origins=*','--user-data-dir='+mkdtempSync(join(tmpdir(),'debelly-review-browser-')),'about:blank'],{stdio:['ignore','ignore','pipe']});
let ws;
try{
 const endpoint=await new Promise((resolve,reject)=>{let log='';chrome.stderr.on('data',data=>{log+=data.toString();const found=log.match(/DevTools listening on (ws:\/\/[^\s]+)/);if(found)resolve(found[1]);});chrome.once('error',reject);setTimeout(()=>reject(new Error('Browser startup timeout')),20000).unref();});
 ws=new WebSocket(endpoint);
 await new Promise((resolve,reject)=>{ws.addEventListener('open',resolve,{once:true});ws.addEventListener('error',reject,{once:true});});
 let next=0; const pending=new Map(), errors=[];
 ws.addEventListener('message',event=>{const d=JSON.parse(event.data);if(d.id){const p=pending.get(d.id);if(p){pending.delete(d.id);if(d.error)p.reject(new Error(JSON.stringify(d.error)));else p.resolve(d.result);}}else if(d.method==='Runtime.exceptionThrown') errors.push(d.params);});
 const send=(method,params={},sessionId)=>new Promise((resolve,reject)=>{const id=++next;pending.set(id,{resolve,reject});ws.send(JSON.stringify({id,method,params,...(sessionId?{sessionId}:{})}));setTimeout(()=>{if(pending.has(id)){pending.delete(id);reject(new Error(method+' timeout'));}},20000).unref();});
 const {targetId}=await send('Target.createTarget',{url:'about:blank'});
 const {sessionId}=await send('Target.attachToTarget',{targetId,flatten:true});
 const call=(method,params)=>send(method,params,sessionId);
 const evaluate=async expression=>{const r=await call('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});if(r.exceptionDetails)throw new Error(JSON.stringify(r.exceptionDetails));return r.result.value;};
 await call('Page.enable');await call('Runtime.enable');
 await call('Emulation.setDeviceMetricsOverride',{width:1280,height:1000,deviceScaleFactor:1,mobile:false});
 await call('Page.navigate',{url:'file://'+resolve(input)});
 await new Promise(r=>setTimeout(r,1200));

 const ready=await evaluate(`(async()=>{const start=Date.now();while(!window.DEBELLY_EVENT_VIEWER||!window.DEBELLY_EVENT_VIEWER.state.samples){if(Date.now()-start>20000)throw new Error('Default segment did not load');await new Promise(r=>setTimeout(r,50));}return window.DEBELLY_EVENT_VIEWER.state;})()`);
 const checks=await evaluate(`(async()=>{const v=window.DEBELLY_EVENT_VIEWER,rows=[];for(const c of CONFIG.cases){for(const segment of c.segments){if(!await v.selectCase(c.key,segment.key))throw new Error('Segment did not load');const p=v.payload;if(v.state.error)throw new Error(v.state.error);v.selectField('all');if(document.querySelectorAll('.plot').length!==7)throw new Error('All-field view missing');for(const key of ['pre_index','post_index','after_step_index']){const button=document.getElementById({pre_index:'beforeButton',post_index:'transitionButton',after_step_index:'afterButton'}[key]);if(p[key]==null){if(!button.disabled)throw new Error('Context event button enabled');}else{button.click();if(v.state.index!==p[key])throw new Error('Event jump mismatch');}}if(segment.event){if(p.samples[p.pre_index].time!==p.event_time||p.samples[p.post_index].time!==p.event_time||p.samples[p.after_step_index].time!==p.event_time+.5)throw new Error('Event stage clock mismatch');}for(const f of CONFIG.fields){v.selectField(f.key);if(document.querySelectorAll('.plot').length!==1)throw new Error('Field selector failed');for(const path of document.querySelectorAll('.current-profile,.reference-profile')){const d=path.getAttribute('d');if(!d||/NaN|Infinity|undefined/.test(d))throw new Error('Invalid SVG profile');}}v.selectField('all');const before=v.state.index;document.getElementById('prevButton').click();document.getElementById('nextButton').click();if(before>0&&v.state.index!==before)throw new Error('Step controls disagree');const rr=document.getElementById('referenceCheck');rr.checked=false;rr.dispatchEvent(new Event('change'));if(document.querySelectorAll('.reference-profile').length)throw new Error('Reference toggle failed');rr.checked=true;rr.dispatchEvent(new Event('change'));const scale=document.getElementById('scaleSelect');scale.value='current';scale.dispatchEvent(new Event('change'));scale.value='window';scale.dispatchEvent(new Event('change'));rows.push({case:c.key,segment:segment.key,samples:p.samples.length,event:p.event_time,all_fields:7});}}await v.selectCase(CONFIG.default_case,CONFIG.default_segment);document.getElementById('afterButton').click();v.selectField('all');return rows;})()`);
 const desktop=await evaluate(`({title:document.title,width:innerWidth,scrollWidth:document.documentElement.scrollWidth,state:window.DEBELLY_EVENT_VIEWER.state,plot_count:document.querySelectorAll('.plot').length,stage:document.getElementById('stageHelp').textContent,links:[...document.querySelectorAll('a')].map(x=>x.getAttribute('href'))})`);
 writeFileSync(join(output,'desktop.png'),Buffer.from((await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false})).data,'base64'));
 await call('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});await evaluate(`window.dispatchEvent(new Event('resize'))`);await new Promise(r=>setTimeout(r,300));
 const mobile=await evaluate(`({width:innerWidth,scrollWidth:document.documentElement.scrollWidth,plot_count:document.querySelectorAll('.plot').length,viewBox:document.querySelector('.plot svg').getAttribute('viewBox')})`);
 writeFileSync(join(output,'mobile.png'),Buffer.from((await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false})).data,'base64'));
 console.log(JSON.stringify({generated_utc:new Date().toISOString(),ready,checks,desktop,mobile,runtime_errors:errors,passed:errors.length===0&&mobile.width===mobile.scrollWidth&&desktop.plot_count===7},null,2));
}finally{if(ws)ws.close();chrome.kill('SIGTERM');}
