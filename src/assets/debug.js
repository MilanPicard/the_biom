// Monitoring des long tasks (blocage UI > 50ms)
try {
    const observer = new PerformanceObserver((list) => {
        for (const entry of list.getEntries()) {
            if (entry.duration > 500) {
                console.warn('🐢 LONG TASK', { duration_ms: entry.duration.toFixed(2) });
            }
        }
    });
    observer.observe({ entryTypes: ['longtask'] });
} catch (e) { }

// Monitoring des appels serveur
const origFetch = window.fetch;
window.fetch = function (url, options) {
    if (url && typeof url === 'string' && url.includes('_dash-update-component')) {
        console.log('📤 CALLBACK REQUEST');
        const start = performance.now();
        return origFetch.apply(this, arguments).then(response => {
            console.log('📥 CALLBACK RESPONSE', { elapsed_ms: (performance.now() - start).toFixed(2) });
            return response;
        });
    }
    return origFetch.apply(this, arguments);
};
