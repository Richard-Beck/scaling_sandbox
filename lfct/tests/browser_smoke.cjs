// NODE_PATH must point to an existing Playwright installation.
const {chromium}=require('playwright');
const path=require('path');
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox']});
 const page=await browser.newPage({viewport:{width:1400,height:1100}});const errors=[];
 page.on('pageerror',e=>errors.push(String(e)));
 await page.goto('file://'+path.resolve(process.argv[2]||'lfct/reports/smoke/index.html'));
 await page.waitForFunction(()=>window.qcTest);
 const images=await page.evaluate(async()=>{const paths=[...new Set(qcTest.data.items.flatMap(i=>i.frames.flatMap(f=>Object.values(f.images))))];for(const p of paths)await new Promise((resolve,reject)=>{const im=new Image();im.onload=()=>im.naturalWidth?resolve():reject(p);im.onerror=()=>reject(p);im.src=p;});return paths.length;});
 await page.click('#next');if(await page.evaluate(()=>qcTest.index)!==1)throw Error('next');
 await page.locator('#slider').evaluate(e=>{e.value='2';e.dispatchEvent(new Event('input'));});
 if(await page.evaluate(()=>qcTest.index)!==2)throw Error('slider');
 await page.click('#plus');if(await page.evaluate(()=>qcTest.zoom)<=1)throw Error('zoom');
 await page.click('#reset');await page.locator('#view').hover();await page.mouse.wheel(0,-35);
 await page.waitForFunction(()=>qcTest.index<2);
 await page.selectOption('#background','rfp');await page.check('#original');await page.check('#unassigned');
 const groups=await page.evaluate(()=>qcTest.data.items.length);
 for(let i=0;i<groups;i++){await page.selectOption('#which',String(i));await page.locator('#view').focus();await page.keyboard.press('End');await page.keyboard.press('Home');}
 await page.selectOption('#which','1');await page.selectOption('#background','phase');
 await page.screenshot({path:process.argv[3]||'lfct/results/smoke/viewer.png',fullPage:true});
 if(errors.length)throw Error(errors.join('\n'));
 console.log(JSON.stringify({status:'passed',images,groups,checks:['file:// assets','next','slider','zoom','wheel','layers','all groups','keyboard'],errors}));await browser.close();
})().catch(e=>{console.error(e);process.exit(1)});
