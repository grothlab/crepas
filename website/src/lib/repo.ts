import { execFileSync } from "node:child_process";
import fs from "node:fs";
import path from "node:path";

// The site is built from inside website/, so the pipeline repo is one level up.
export const REPO_ROOT = process.env.CREPAS_REPO_ROOT ?? path.resolve(process.cwd(), "..");

// Public address of the site; absolute links to it in the docs are rewritten to the version being viewed.
export const SITE_URL = "https://crepas.grothlab.workers.dev";

function git(args: string[]): string {
    return execFileSync("git", args, { cwd: REPO_ROOT, stdio: ["ignore", "pipe", "ignore"] }).toString();
}

// `ref` null reads the checkout (the development version); a tag reads that release.
export function readRepoFile(relPath: string, ref: string | null = null): string {
    if (ref === null) return fs.readFileSync(path.join(REPO_ROOT, relPath), "utf8");
    return git(["show", `${ref}:${relPath}`]);
}

export function repoFileExists(relPath: string, ref: string | null = null): boolean {
    if (ref === null) return fs.existsSync(path.join(REPO_ROOT, relPath));
    try {
        git(["cat-file", "-e", `${ref}:${relPath}`]);
        return true;
    } catch {
        return false;
    }
}

export interface Contributor {
    name: string;
    affiliation?: string;
    github?: string;
    orcid?: string;
    contribution: string[];
}

export interface Manifest {
    name: string;
    description: string;
    homePage: string;
    version: string;
    nextflowVersion: string;
    contributors: Contributor[];
}

function quotedValue(text: string, key: string, separator: "=" | ":"): string | undefined {
    const re = new RegExp(`\\b${key}\\s*${separator}\\s*(?:"""([\\s\\S]*?)"""|'([^']*)'|"([^"]*)")`);
    const match = text.match(re);
    return match?.slice(1).find((v) => v !== undefined)?.trim();
}

const manifestCache = new Map<string, Manifest>();

export function getManifest(ref: string | null = null): Manifest {
    const key = ref ?? "";
    const cached = manifestCache.get(key);
    if (cached) return cached;
    const config = readRepoFile("nextflow.config", ref);
    const block = config.match(/^manifest\s*\{([\s\S]*?)^\}/m)?.[1] ?? "";
    const contributorsText = block.match(/contributors\s*=\s*\[([\s\S]*?)\n\s*\]\s*\n/)?.[1] ?? "";
    const contributors = contributorsText
        .split(/\bname\s*:/)
        .slice(1)
        .map((chunk) => {
            const text = "name:" + chunk;
            const contribution = text.match(/contribution\s*:\s*\[([^\]]*)\]/)?.[1] ?? "";
            return {
                name: quotedValue(text, "name", ":") ?? "",
                affiliation: quotedValue(text, "affiliation", ":"),
                github: quotedValue(text, "github", ":")?.replace(/^@/, ""),
                orcid: quotedValue(text, "orcid", ":"),
                contribution: [...contribution.matchAll(/'([^']+)'/g)].map((m) => m[1]),
            };
        })
        .filter((c) => c.name);

    const manifest = {
        name: quotedValue(block, "name", "=") ?? "grothlab/crepas",
        description: quotedValue(block, "description", "=") ?? "",
        homePage: quotedValue(block, "homePage", "=") ?? "https://github.com/grothlab/crepas",
        version: quotedValue(block, "version", "=") ?? "dev",
        nextflowVersion: quotedValue(block, "nextflowVersion", "=") ?? "",
        contributors,
    };
    manifestCache.set(key, manifest);
    return manifest;
}

// Cloudflare Workers Builds exposes the branch being built; fall back to git locally.
export function getBranch(): string {
    if (process.env.WORKERS_CI_BRANCH) return process.env.WORKERS_CI_BRANCH;
    try {
        const branch = git(["rev-parse", "--abbrev-ref", "HEAD"]).trim();
        if (branch && branch !== "HEAD") return branch;
    } catch {
        // not a git checkout
    }
    return "dev";
}

export interface SiteVersion {
    // URL segment and dropdown label: a release tag, or "dev"
    id: string;
    // git ref to read files from; null for the checkout
    ref: string | null;
    isDev: boolean;
    isLatest: boolean;
    // tag date (YYYY-MM-DD) for releases
    date?: string;
}

function compareSemverDesc(a: string, b: string): number {
    const pa = a.split(".").map(Number);
    const pb = b.split(".").map(Number);
    for (let i = 0; i < 3; i++) if (pa[i] !== pb[i]) return pb[i] - pa[i];
    return 0;
}

let versionsCache: SiteVersion[] | undefined;

// Releases are the X.Y.Z tags in the checkout (run `git fetch --tags` before building), newest first, then dev.
export function getVersions(): SiteVersion[] {
    if (versionsCache) return versionsCache;
    let tags: string[] = [];
    try {
        tags = git(["tag", "--list"])
            .split("\n")
            .map((t) => t.trim())
            .filter((t) => /^\d+\.\d+\.\d+$/.test(t))
            .sort(compareSemverDesc);
    } catch {
        // not a git checkout: only the development version is built
    }
    const releases: SiteVersion[] = tags.map((tag, i) => ({
        id: tag,
        ref: tag,
        isDev: false,
        isLatest: i === 0,
        date: git(["log", "-1", "--format=%cs", tag]).trim(),
    }));
    versionsCache = [...releases, { id: "dev", ref: null, isDev: true, isLatest: releases.length === 0 }];
    return versionsCache;
}

export function getLatestVersion(): SiteVersion {
    return getVersions().find((v) => v.isLatest)!;
}

export function getVersion(id: string): SiteVersion {
    return getVersions().find((v) => v.id === id)!;
}

// Ref used in GitHub links for a version: the tag for releases, the checked-out branch for dev.
export function githubRef(version: SiteVersion): string {
    return version.ref ?? getBranch();
}

// Date of the last commit on the checkout, shown as "last update".
export function getLastUpdate(): string | undefined {
    try {
        return git(["log", "-1", "--format=%cs"]).trim() || undefined;
    } catch {
        return undefined;
    }
}

export const PAGES = ["", "usage", "parameters", "output", "releases"] as const;
export type Page = (typeof PAGES)[number];

export function pageHref(version: SiteVersion, page: Page): string {
    return `/${version.id}/${page ? `${page}/` : ""}`;
}

export interface Components {
    modules: string[];
    subworkflows: string[];
}

// nf-core modules and subworkflows installed in a version, from its modules.json.
export function getComponents(ref: string | null): Components {
    if (!repoFileExists("modules.json", ref)) return { modules: [], subworkflows: [] };
    const json = JSON.parse(readRepoFile("modules.json", ref));
    const modules = new Set<string>();
    const subworkflows = new Set<string>();
    for (const repo of Object.values(json.repos ?? {}) as Record<string, Record<string, object>>[]) {
        for (const org of Object.values(repo.modules ?? {})) Object.keys(org).forEach((m) => modules.add(m));
        for (const org of Object.values(repo.subworkflows ?? {})) Object.keys(org).forEach((s) => subworkflows.add(s));
    }
    return { modules: [...modules].sort(), subworkflows: [...subworkflows].sort() };
}
