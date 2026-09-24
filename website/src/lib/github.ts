// Repository statistics from the GitHub API, fetched once at build time. A failed or
// rate-limited request leaves the statistics out rather than failing the build.
// Set GITHUB_TOKEN to raise the API rate limit.

export interface RepoStats {
    subscribers: number;
    stars: number;
    forks: number;
    openIssues: number;
    openPulls: number;
    contributors?: number;
}

async function api(url: string): Promise<Response> {
    const headers: Record<string, string> = { Accept: "application/vnd.github+json" };
    if (process.env.GITHUB_TOKEN) headers.Authorization = `Bearer ${process.env.GITHUB_TOKEN}`;
    const response = await fetch(url, { headers, signal: AbortSignal.timeout(10000) });
    if (!response.ok) throw new Error(`${url}: HTTP ${response.status}`);
    return response;
}

// Number of items in a paginated listing, read from the last-page link with one item per page.
async function countItems(url: string): Promise<number> {
    const response = await api(`${url}${url.includes("?") ? "&" : "?"}per_page=1`);
    const last = response.headers.get("link")?.match(/[?&]page=(\d+)[^>]*>;\s*rel="last"/);
    if (last) return Number(last[1]);
    return ((await response.json()) as unknown[]).length;
}

let statsPromise: Promise<RepoStats | undefined> | undefined;

export function getRepoStats(homePage: string): Promise<RepoStats | undefined> {
    statsPromise ??= (async () => {
        const slug = homePage.replace(/^https:\/\/github\.com\//, "").replace(/\/$/, "");
        const base = `https://api.github.com/repos/${slug}`;
        try {
            const repo = (await (await api(base)).json()) as {
                subscribers_count: number;
                stargazers_count: number;
                forks_count: number;
                open_issues_count: number;
            };
            // open_issues_count includes pull requests
            const openPulls = await countItems(`${base}/pulls?state=open`);
            const contributors = await countItems(`${base}/contributors?anon=1`).catch(() => undefined);
            return {
                subscribers: repo.subscribers_count,
                stars: repo.stargazers_count,
                forks: repo.forks_count,
                openIssues: repo.open_issues_count - openPulls,
                openPulls,
                contributors,
            };
        } catch (error) {
            console.warn(`[github] repository statistics left out: ${(error as Error).message}`);
            return undefined;
        }
    })();
    return statsPromise;
}
