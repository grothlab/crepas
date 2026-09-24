import { renderInline, renderMarkdown } from "./markdown";
import { readRepoFile, type SiteVersion } from "./repo";

interface RawProperty {
    type?: string | string[];
    description?: string;
    help_text?: string;
    fa_icon?: string;
    hidden?: boolean;
    default?: unknown;
    enum?: unknown[];
    pattern?: string;
    format?: string;
    minimum?: number;
    maximum?: number;
}

interface RawGroup {
    title?: string;
    description?: string;
    help_text?: string;
    fa_icon?: string;
    required?: string[];
    properties?: Record<string, RawProperty>;
}

export interface SchemaProperty {
    name: string;
    type: string;
    faIcon?: string;
    hidden: boolean;
    required: boolean;
    descriptionHtml: string;
    helpTextHtml: string;
    default?: string;
    enumValues?: string[];
    pattern?: string;
    format?: string;
    minimum?: number;
    maximum?: number;
}

export interface SchemaGroup {
    id: string;
    title: string;
    faIcon?: string;
    descriptionHtml: string;
    helpTextHtml: string;
    hidden: boolean;
    properties: SchemaProperty[];
}

async function toProperty(name: string, raw: RawProperty, required: string[], version: SiteVersion): Promise<SchemaProperty> {
    return {
        name,
        type: Array.isArray(raw.type) ? raw.type.join(" | ") : (raw.type ?? ""),
        faIcon: raw.fa_icon,
        hidden: raw.hidden === true,
        required: required.includes(name),
        descriptionHtml: await renderInline(raw.description, version),
        helpTextHtml: raw.help_text
            ? (await renderMarkdown(raw.help_text, { sourceFile: "nextflow_schema.json", version, highlight: false })).html
            : "",
        default: raw.default === undefined || raw.default === "" ? undefined : String(raw.default),
        enumValues: raw.enum?.map(String),
        pattern: raw.pattern,
        format: raw.format,
        minimum: raw.minimum,
        maximum: raw.maximum,
    };
}

export async function getSchemaGroups(version: SiteVersion): Promise<SchemaGroup[]> {
    const schema = JSON.parse(readRepoFile("nextflow_schema.json", version.ref));
    const defs: Record<string, RawGroup> = schema.$defs ?? schema.definitions ?? {};
    const ordered: string[] = (schema.allOf ?? [])
        .map((entry: { $ref?: string }) => entry.$ref?.split("/").pop())
        .filter((id: string | undefined): id is string => Boolean(id && defs[id]));
    const ids = [...ordered, ...Object.keys(defs).filter((id) => !ordered.includes(id))];

    const groups = await Promise.all(
        ids.map(async (id) => {
            const group = defs[id];
            const required = group.required ?? [];
            const properties = await Promise.all(
                Object.entries(group.properties ?? {}).map(([name, raw]) => toProperty(name, raw, required, version)),
            );
            return {
                id: id.replaceAll("_", "-"),
                title: group.title ?? id,
                faIcon: group.fa_icon,
                descriptionHtml: await renderInline(group.description, version),
                helpTextHtml: group.help_text
                    ? (await renderMarkdown(group.help_text, { sourceFile: "nextflow_schema.json", version, highlight: false })).html
                    : "",
                hidden: properties.length > 0 && properties.every((p) => p.hidden),
                properties,
            };
        }),
    );

    // Parameters declared outside any group still deserve a listing.
    if (schema.properties && Object.keys(schema.properties).length > 0) {
        const required: string[] = schema.required ?? [];
        const properties = await Promise.all(
            Object.entries(schema.properties as Record<string, RawProperty>).map(([name, raw]) =>
                toProperty(name, raw, required, version),
            ),
        );
        groups.push({
            id: "other-parameters",
            title: "Other parameters",
            faIcon: "fas fa-ellipsis",
            descriptionHtml: "",
            helpTextHtml: "",
            hidden: properties.every((p) => p.hidden),
            properties,
        });
    }
    return groups;
}
