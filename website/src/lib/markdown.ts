import path from "node:path";
import rehypeShiki from "@shikijs/rehype";
import type { Element, Root as HastRoot } from "hast";
import { toString } from "hast-util-to-string";
import type { Blockquote, Paragraph, Root as MdastRoot } from "mdast";
import rehypeRaw from "rehype-raw";
import rehypeSlug from "rehype-slug";
import rehypeStringify from "rehype-stringify";
import remarkGfm from "remark-gfm";
import remarkParse from "remark-parse";
import remarkRehype from "remark-rehype";
import { unified } from "unified";
import { visit } from "unist-util-visit";
import { getManifest, githubRef, type Page, pageHref, readRepoFile, SITE_URL, type SiteVersion } from "./repo";

export interface Heading {
    depth: number;
    slug: string;
    text: string;
}

// Repo files that have their own page on the site; anything else links to GitHub.
const PAGE_FOR_FILE: Record<string, Page> = {
    "README.md": "",
    "docs/usage.md": "usage",
    "docs/output.md": "output",
    "CHANGELOG.md": "releases",
    "nextflow_schema.json": "parameters",
};

const ALERTS: Record<string, { variant: string; icon: string; title: string }> = {
    NOTE: { variant: "info", icon: "fa-circle-info", title: "Note" },
    TIP: { variant: "success", icon: "fa-lightbulb", title: "Tip" },
    IMPORTANT: { variant: "primary", icon: "fa-circle-exclamation", title: "Important" },
    WARNING: { variant: "warning", icon: "fa-triangle-exclamation", title: "Warning" },
    CAUTION: { variant: "danger", icon: "fa-circle-xmark", title: "Caution" },
};

// GitHub `> [!NOTE]` blockquotes become Bootstrap alerts, as on nf-co.re.
function remarkGithubAlerts() {
    return (tree: MdastRoot) => {
        visit(tree, "blockquote", (node: Blockquote) => {
            const first = node.children[0];
            if (first?.type !== "paragraph") return;
            const text = first.children[0];
            if (text?.type !== "text") return;
            const match = text.value.match(/^\[!(NOTE|TIP|IMPORTANT|WARNING|CAUTION)\]\s*/);
            if (!match) return;

            text.value = text.value.slice(match[0].length);
            if (!text.value) first.children.shift();
            if (first.children[0]?.type === "break") first.children.shift();
            if (first.children.length === 0) node.children.shift();

            const alert = ALERTS[match[1]];
            node.data = {
                hName: "div",
                hProperties: { className: ["alert", `alert-${alert.variant}`, "markdown-alert"], role: "alert" },
            };
            const heading: Paragraph = {
                type: "paragraph",
                data: { hProperties: { className: ["alert-heading", "fw-semibold", "mb-2"] } },
                children: [
                    { type: "html", value: `<i class="fa-solid ${alert.icon} me-2"></i>` },
                    { type: "text", value: alert.title },
                ],
            };
            node.children.unshift(heading);
        });
    };
}

function mdastToString(node: { value?: unknown; children?: unknown[] }): string {
    if (typeof node.value === "string") return node.value;
    return (node.children ?? []).map((child) => mdastToString(child as typeof node)).join("");
}

// The page header already shows the pipeline name and badges add noise, so drop both.
// Badges only sit above the first section heading; image-only paragraphs further down
// (e.g. the metro map) are content and stay.
function remarkStripTitleAndBadges() {
    return (tree: MdastRoot) => {
        const firstH1 = tree.children.findIndex((n) => n.type === "heading" && n.depth === 1);
        if (firstH1 !== -1) tree.children.splice(firstH1, 1);
        const firstSection = tree.children.findIndex((n) => n.type === "heading");
        tree.children = tree.children.filter((node, index) => {
            if (firstSection !== -1 && index >= firstSection) return true;
            // The docs' "Please read this documentation on …" pointer is redundant on the site itself.
            if (node.type === "blockquote") return !/Please read this documentation/.test(mdastToString(node));
            if (node.type !== "paragraph") return true;
            const onlyBadges = node.children.every(
                (child) =>
                    child.type === "image" ||
                    (child.type === "link" && child.children.every((c) => c.type === "image")) ||
                    (child.type === "text" && !child.value.trim()),
            );
            return !onlyBadges;
        });
    };
}

// Keeps links within the version being viewed: repo files with a page, and absolute links to the
// site, go to that page of the same version; other repo files link to GitHub at the version's ref.
function rehypeRepoLinks(options: { sourceFile: string; version: SiteVersion }) {
    const { version } = options;
    const { homePage } = getManifest(version.ref);
    const ref = githubRef(version);
    const sourceDir = path.posix.dirname(options.sourceFile);
    return (tree: HastRoot) => {
        visit(tree, "element", (node: Element) => {
            const attr = node.tagName === "a" ? "href" : node.tagName === "img" ? "src" : undefined;
            if (!attr) return;
            const url = node.properties?.[attr];
            if (typeof url !== "string" || !url) return;

            if (attr === "href" && url.startsWith(`${SITE_URL}/`)) {
                const [sitePath, hash] = url.slice(SITE_URL.length + 1).split("#");
                const page = sitePath.replace(/\/$/, "") as Page;
                if ((Object.values(PAGE_FOR_FILE) as string[]).includes(page)) {
                    node.properties[attr] = pageHref(version, page) + (hash ? `#${hash}` : "");
                }
                return;
            }
            if (/^([a-z][a-z0-9+.-]*:|#|\/)/i.test(url)) return;

            const [target, hash] = url.split("#");
            const resolved = path.posix.normalize(path.posix.join(sourceDir, target));
            const page = PAGE_FOR_FILE[resolved];
            if (attr === "href" && page !== undefined) {
                node.properties[attr] = pageHref(version, page) + (hash ? `#${hash}` : "");
                return;
            }
            const suffix = attr === "src" ? "?raw=true" : hash ? `#${hash}` : "";
            node.properties[attr] = `${homePage}/blob/${ref}/${resolved}${suffix}`;
        });
    };
}

function rehypeBootstrap(options: { headings?: Heading[] }) {
    return (tree: HastRoot) => {
        visit(tree, "element", (node: Element, index, parent) => {
            if (node.tagName === "table") {
                node.properties.className = ["table", "table-sm", "table-hover"];
                // Wide tables (e.g. the example samplesheets) scroll inside the column instead of overflowing it.
                if (parent && index !== undefined && !(parent.type === "element" && parent.tagName === "div")) {
                    parent.children[index] = {
                        type: "element",
                        tagName: "div",
                        properties: { className: ["table-responsive"] },
                        children: [node],
                    };
                }
            }
            if (options.headings && (node.tagName === "h2" || node.tagName === "h3")) {
                const slug = node.properties?.id;
                if (typeof slug === "string") {
                    options.headings.push({ depth: Number(node.tagName[1]), slug, text: toString(node) });
                }
            }
        });
    };
}

export async function renderMarkdown(
    source: string,
    options: { sourceFile: string; version: SiteVersion; stripTitle?: boolean; highlight?: boolean },
): Promise<{ html: string; headings: Heading[] }> {
    const headings: Heading[] = [];
    const processor = unified().use(remarkParse).use(remarkGfm).use(remarkGithubAlerts);
    if (options.stripTitle) processor.use(remarkStripTitleAndBadges);
    processor
        .use(remarkRehype, { allowDangerousHtml: true })
        .use(rehypeRaw)
        .use(rehypeSlug)
        .use(rehypeRepoLinks, { sourceFile: options.sourceFile, version: options.version })
        .use(rehypeBootstrap, { headings });
    if (options.highlight !== false) {
        processor.use(rehypeShiki, {
            themes: { light: "github-light", dark: "github-dark" },
            defaultColor: false,
            fallbackLanguage: "text",
            lazy: true,
        });
    }
    const file = await processor.use(rehypeStringify).process(source);
    return { html: String(file), headings };
}

export async function renderRepoMarkdown(relPath: string, version: SiteVersion, stripTitle = false) {
    return renderMarkdown(readRepoFile(relPath, version.ref), { sourceFile: relPath, version, stripTitle });
}

// Short schema descriptions render inline, without the wrapping <p>.
export async function renderInline(source: string | undefined, version: SiteVersion): Promise<string> {
    if (!source) return "";
    const { html } = await renderMarkdown(source, { sourceFile: "nextflow_schema.json", version, highlight: false });
    return html.replace(/^<p>([\s\S]*)<\/p>\s*$/, "$1");
}
