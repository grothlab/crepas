import { readRepoFile } from "./repo";

export interface ChangelogRelease {
    version: string;
    name?: string;
    date?: string;
    // changelog section without its heading
    body: string;
}

// Released versions from the checkout's CHANGELOG.md (it holds every release), newest first.
// Headings look like `## [[1.0.0](https://…/releases/tag/1.0.0)] - Mercurian Cinnabar - 2026-06-21`;
// development sections (`1.1.0dev`) are left out.
export function getChangelogReleases(): ChangelogRelease[] {
    const changelog = readRepoFile("CHANGELOG.md");
    const sections = changelog.split(/^(?=## )/m).filter((s) => s.startsWith("## "));
    return sections
        .map((section) => {
            const [heading, ...rest] = section.split("\n");
            const version = heading.match(/\d+\.\d+\.\d+[0-9A-Za-z.-]*/)?.[0] ?? "";
            const parts = heading.replace(/^##\s*/, "").split(/\s+-\s+/);
            const date = heading.match(/\d{4}-\d{2}-\d{2}/)?.[0];
            const name = parts.length >= 3 ? parts[1].trim() : undefined;
            return { version, name, date, body: rest.join("\n").trim() };
        })
        .filter((release) => /^\d+\.\d+\.\d+$/.test(release.version));
}
